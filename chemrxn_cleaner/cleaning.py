# chemrxn_cleaner/cleaning.py

from __future__ import annotations

from typing import Iterable, List, Optional, Tuple, Dict, Any

from .types import ReactionRecord
from .parsing import parse_reaction_smiles, canonicalize_reaction
from .filters import ReactionFilter, default_filters


def clean_reactions(
    rxn_smiles_list: Iterable[Tuple[str, Dict[str, Any]]],
    filters: Optional[List[ReactionFilter]] = None,
    drop_failed_parse: bool = True,
    strict: bool = True,
) -> List[ReactionRecord]:
    """
    Parse and clean a list of reaction SMILES.

    This is the core entry point for the cleaning pipeline.

    Args:
        rxn_smiles_list:
            An iterable of (reaction_smiles, meta_dict) tuples.

        filters:
            A list of predicate functions. Each filter takes a ReactionRecord and
            returns True if the reaction should be kept, False otherwise.
            If None, uses `default_filters()`.

        drop_failed_parse:
            - If True: silently drop reactions that cannot be parsed.
            - If False: re-raise the exception from parsing.

        strict:
            Passed to `parse_reaction_smiles`:
            - True: require exactly 3 parts ('reactants>reagents>products').
            - False: auto-pad/truncate to 3 parts.

    Returns:
        A list of cleaned ReactionRecord objects which passed all filters.
    """
    if filters is None:
        filters = default_filters()

    cleaned: List[ReactionRecord] = []

    for rxn_entry in rxn_smiles_list:
        if rxn_entry is None:
            continue
        rxn, meta = rxn_entry
        if rxn is None:
            continue
        if meta is None:
            meta = {}

        try:
            record = parse_reaction_smiles(rxn, strict=strict)
        except Exception:
            if drop_failed_parse:
                continue
            else:
                raise
        record.meta = meta

        keep = True
        for f in filters:
            if not f(record):
                keep = False
                break

        if keep:
            cleaned.append(record)

    return cleaned


def clean_and_canonicalize(
    rxn_smiles_list: Iterable[Tuple[str, Dict[str, Any]]],
    filters: Optional[List[ReactionFilter]] = None,
    drop_failed_parse: bool = True,
    strict: bool = True,
    isomeric: bool = True,
) -> List[ReactionRecord]:
    """
    Clean reactions and canonicalize all SMILES in one pass.

    This is a convenience wrapper around `clean_reactions` + `canonicalize_reaction`.

    Args:
        rxn_smiles_list:
            Iterable of (reaction_smiles, meta_dict) pairs.

        filters:
            List of ReactionFilter predicates. If None, uses default_filters().

        drop_failed_parse:
            Whether to silently drop unparseable reactions.

        strict:
            Whether to enforce exactly three '>' parts in reaction SMILES.

        isomeric:
            Whether to keep isomeric SMILES when canonicalizing.

    Returns:
        List[ReactionRecord] with canonicalized reactants/reagents/products.
    """
    cleaned = clean_reactions(
        rxn_smiles_list=rxn_smiles_list,
        filters=filters,
        drop_failed_parse=drop_failed_parse,
        strict=strict,
    )

    canon_records: List[ReactionRecord] = []
    for rec in cleaned:
        canon_records.append(canonicalize_reaction(rec, isomeric=isomeric))

    return canon_records


def basic_cleaning_pipeline(
    rxn_smiles_list: Iterable[Tuple[str, Dict[str, Any]]],
) -> List[ReactionRecord]:
    """
    A simple out-of-the-box cleaning pipeline.

    Equivalent to:
        clean_and_canonicalize(
            rxn_smiles_list,
            filters=default_filters(),
            drop_failed_parse=True,
            strict=True,
            isomeric=True,
        )
    """
    return clean_and_canonicalize(
        rxn_smiles_list=rxn_smiles_list,
        filters=default_filters(),
        drop_failed_parse=True,
        strict=True,
        isomeric=True,
    )
