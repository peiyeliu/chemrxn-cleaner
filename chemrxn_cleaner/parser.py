# chemrxn_cleaner/parser.py

from __future__ import annotations

import logging
import re
from typing import Iterable, List, Optional, Tuple

from rdkit import Chem

from .types import ReactionRecord

logger = logging.getLogger(__name__)

# Regex to remove atom-map numbers inside bracketed atoms, e.g., ``[CH3:1]``.
_ATOM_MAP_PATTERN = re.compile(r"(\[[^\[\]]*?):\d+([^\[\]]*\])")

# ---------------------- basic SMILES tool functions ---------------------- #


def parse_smiles_list(smiles_block: str) -> List[str]:
    """Split a dot-delimited SMILES block into individual tokens.

    Args:
        smiles_block: String such as ``'CCO.CCBr'``.

    Returns:
        List of trimmed SMILES strings; empty input yields an empty list.
    """
    if not smiles_block:
        return []

    parts = [s.strip() for s in smiles_block.split(".")]
    return [s for s in parts if s]


def canonicalize_smiles(smiles: str, isomeric: bool = True) -> str:
    """Convert a SMILES string to its canonical form using RDKit.

    Args:
        smiles: Input SMILES string.
        isomeric: Preserve stereochemistry in the canonical output.

    Returns:
        Canonical SMILES string.

    Raises:
        ValueError: If the SMILES string is empty or cannot be parsed.
    """
    if not smiles:
        logger.error("Empty SMILES string cannot be canonicalized.")
        raise ValueError("Empty SMILES string cannot be canonicalized.")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.error("Invalid SMILES (cannot parse): %r", smiles)
        raise ValueError(f"Invalid SMILES (cannot parse): {smiles!r}")

    return Chem.MolToSmiles(mol, isomericSmiles=isomeric)


def canonicalize_smiles_list(
    smiles_list: Iterable[str],
    isomeric: bool = True,
) -> List[str]:
    """Canonicalize each SMILES string in a list.

    Args:
        smiles_list: Iterable of SMILES strings to convert.
        isomeric: Whether to preserve stereochemistry in each output SMILES.

    Returns:
        List of canonical SMILES strings.

    Raises:
        ValueError: If any SMILES cannot be parsed.
    """
    return [canonicalize_smiles(s, isomeric=isomeric) for s in smiles_list]


# ---------------------- Reaction SMILES parser ---------------------- #


def _strip_atom_mapping_tokens(rxn_smiles: str) -> Tuple[str, Optional[str]]:
    """Remove atom-map numbers from a reaction SMILES string.

    Args:
        rxn_smiles: Reaction SMILES that may contain atom-map numbers such as
            ``[CH3:1]``.

    Returns:
        Tuple of ``(stripped_smiles, original_if_stripped)`` where the second
        element is ``None`` when no atom-maps were removed.
    """

    original = rxn_smiles
    updated = rxn_smiles
    while True:
        stripped = _ATOM_MAP_PATTERN.sub(r"\1\2", updated)
        if stripped == updated:
            break
        updated = stripped

    if updated == original:
        return updated, None
    return updated, original


def parse_reaction_smiles(
    rxn_smiles: str,
    strict: bool = True,
    strip_atom_mapping: bool = False,
) -> ReactionRecord:
    """Parse a reaction SMILES string into a ``ReactionRecord``.

    Reaction SMILES must follow ``reactants>reagents>products``. When
    ``strict`` is False, missing fields are padded with empty strings instead
    of raising.

    Args:
        rxn_smiles: Raw reaction SMILES string to parse.
        strict: Enforce exactly three ``>``-separated parts when True.
        strip_atom_mapping: When True, remove atom-map numbers and store the
            original mapped string in ``atom_mapping``.

    Returns:
        ReactionRecord with parsed reactants, reagents, and products lists.

    Raises:
        ValueError: If ``rxn_smiles`` is None or incorrectly formatted when
            ``strict`` is True.
    """
    if rxn_smiles is None:
        logger.error("rxn_smiles cannot be None.")
        raise ValueError("rxn_smiles cannot be None.")

    rxn_smiles = rxn_smiles.strip()

    atom_mapping: Optional[str] = None
    if strip_atom_mapping:
        rxn_smiles, atom_mapping = _strip_atom_mapping_tokens(rxn_smiles)

    parts = rxn_smiles.split(">")
    if len(parts) != 3:
        if strict:
            logger.error("Invalid reaction SMILES format: %r", rxn_smiles)
            raise ValueError(f"Invalid reaction SMILES format: {rxn_smiles!r}")
        else:
            # Pad or truncate to length 3
            parts = (parts + ["", "", ""])[:3]

    reactants_block, reagents_block, products_block = parts

    reactants = parse_smiles_list(reactants_block)
    reagents = parse_smiles_list(reagents_block)
    products = parse_smiles_list(products_block)

    return ReactionRecord(
        reaction_smiles=rxn_smiles,
        reactants=reactants,
        reagents=reagents,
        products=products,
        atom_mapping=atom_mapping,
    )


def canonicalize_reaction(
    record: ReactionRecord,
    isomeric: bool = True,
) -> ReactionRecord:
    """Return a new ReactionRecord with all SMILES canonicalized.

    Args:
        record: Reaction record to convert.
        isomeric: Preserve stereochemistry in canonical SMILES when True.

    Returns:
        A new ReactionRecord containing canonicalized SMILES and copied metadata.
    """
    canon_reactants = canonicalize_smiles_list(record.reactants, isomeric=isomeric)
    canon_reagents = canonicalize_smiles_list(record.reagents, isomeric=isomeric)
    canon_products = canonicalize_smiles_list(record.products, isomeric=isomeric)

    return ReactionRecord(
        reaction_smiles=record.reaction_smiles,
        reactants=canon_reactants,
        reagents=canon_reagents,
        products=canon_products,
        extra_metadata=dict(record.extra_metadata or {}),
    )
