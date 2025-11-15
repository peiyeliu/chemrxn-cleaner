# chemrxn_cleaner/parsing.py

from __future__ import annotations

from typing import Iterable, List

from rdkit import Chem

from .types import ReactionRecord


# ---------------------- basic SMILES tool functions ---------------------- #


def parse_smiles_list(smiles_block: str) -> List[str]:
    """
    Split a SMILES block like 'CCO.CCBr' into ['CCO', 'CCBr'].

    - Empty string -> []
    - Strips whitespace around each token
    - Skips empty tokens
    """
    if not smiles_block:
        return []

    parts = [s.strip() for s in smiles_block.split(".")]
    return [s for s in parts if s]


def canonicalize_smiles(smiles: str, isomeric: bool = True) -> str:
    """
    Convert a SMILES string to its canonical form using RDKit.

    Raises:
        ValueError if SMILES cannot be parsed.
    """
    if not smiles:
        raise ValueError("Empty SMILES string cannot be canonicalized.")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES (cannot parse): {smiles!r}")

    return Chem.MolToSmiles(mol, isomericSmiles=isomeric)


def canonicalize_smiles_list(
    smiles_list: Iterable[str],
    isomeric: bool = True,
) -> List[str]:
    """
    Canonicalize a list of SMILES strings.

    Invalid SMILES will raise ValueError.
    """
    return [canonicalize_smiles(s, isomeric=isomeric) for s in smiles_list]


# ---------------------- Reaction SMILES parser ---------------------- #


def parse_reaction_smiles(
    rxn_smiles: str,
    strict: bool = True,
) -> ReactionRecord:
    """
    Parse a reaction SMILES of the form:

        reactants>reagents>products

    Examples:
        "C=CCBr>>C=CCI"
        ""CC(=O)O.OCC>[H+].[Cl-].OCC>CC(=O)OCC"
        "CC(=O)Cl.NH3>>CC(=O)NH2"

    Args:
        rxn_smiles: Raw reaction SMILES string.
        strict:
            - If True  (default): require exactly three '>'-separated parts,
              otherwise raise ValueError.
            - If False: if parts are fewer than 3, pad missing fields with ''.

    Returns:
        ReactionRecord with reactants/reagents/products as lists of SMILES.

    Raises:
        ValueError on invalid format when strict=True.
    """
    if rxn_smiles is None:
        raise ValueError("rxn_smiles cannot be None.")

    parts = rxn_smiles.strip().split(">")
    if len(parts) != 3:
        if strict:
            raise ValueError(f"Invalid reaction SMILES format: {rxn_smiles!r}")
        else:
            # Pad or truncate to length 3
            parts = (parts + ["", "", ""])[:3]

    reactants_block, reagents_block, products_block = parts

    reactants = parse_smiles_list(reactants_block)
    reagents = parse_smiles_list(reagents_block)
    products = parse_smiles_list(products_block)

    return ReactionRecord(
        raw=rxn_smiles,
        reactants=reactants,
        reagents=reagents,
        products=products,
        meta=None,
    )


def canonicalize_reaction(
    record: ReactionRecord,
    isomeric: bool = True,
) -> ReactionRecord:
    """
    Return a new ReactionRecord with all SMILES canonicalized.

    This does NOT modify the input record in-place.
    """
    canon_reactants = canonicalize_smiles_list(record.reactants, isomeric=isomeric)
    canon_reagents = canonicalize_smiles_list(record.reagents, isomeric=isomeric)
    canon_products = canonicalize_smiles_list(record.products, isomeric=isomeric)

    return ReactionRecord(
        raw=record.raw,
        reactants=canon_reactants,
        reagents=canon_reagents,
        products=canon_products,
        meta=record.meta,
    )
