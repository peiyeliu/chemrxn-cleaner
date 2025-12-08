# chemrxn_cleaner/cleaner.py

from __future__ import annotations

import logging
from typing import Iterable, List, Optional, Tuple

from .filters import ReactionFilter, default_filters
from .parser import canonicalize_reaction, parse_reaction_smiles
from .reporter import CleaningStats, FilterStats
from .types import ReactionRecord

logger = logging.getLogger(__name__)


def _get_filter_name(f: ReactionFilter) -> str:
    return getattr(f, "__name__", f.__class__.__name__)


def _ensure_filter_stats(stats: CleaningStats, filter_name: str) -> FilterStats:
    if filter_name not in stats.per_filter:
        stats.per_filter[filter_name] = FilterStats(name=filter_name)
    return stats.per_filter[filter_name]


def _clean_reactions_internal(
    rxn_smiles_list: Iterable[ReactionRecord],
    filters: Optional[List[ReactionFilter]],
    drop_failed_parse: bool,
    strict: bool,
    collect_stats: bool,
) -> Tuple[List[ReactionRecord], CleaningStats]:
    if filters is None:
        filters = default_filters()

    logger.info(
        "Starting cleaning pipeline (drop_failed_parse=%s, strict=%s, filters=%d)",
        drop_failed_parse,
        strict,
        len(filters),
    )

    cleaned: List[ReactionRecord] = []
    stats = CleaningStats()

    for idx, rxn_entry in enumerate(rxn_smiles_list):
        if rxn_entry is None:
            logger.debug("Skipping None reaction entry at index %d", idx)
            continue

        if collect_stats:
            stats.n_input += 1

        record = rxn_entry
        if record.extra_metadata is None:
            record.extra_metadata = {}

        needs_parse = not (record.reactants or record.reagents or record.products)
        if needs_parse:
            try:
                parsed = parse_reaction_smiles(record.reaction_smiles, strict=strict)
            except Exception as exc:
                if collect_stats:
                    stats.n_failed_parse += 1
                if drop_failed_parse:
                    logger.warning(
                        "Dropping reaction at index %d: failed to parse reaction"
                        "(%s)",
                        idx,
                        exc,
                    )
                    continue
                raise
            record.reaction_smiles = parsed.reaction_smiles
            record.reactants = parsed.reactants
            record.reagents = parsed.reagents
            record.products = parsed.products

        keep = True
        failing_filter = None
        for f in filters:
            filter_name = _get_filter_name(f)
            if collect_stats:
                fstats = _ensure_filter_stats(stats, filter_name)
                fstats.applied += 1
            passed = f(record)
            if passed:
                if collect_stats:
                    fstats.passed += 1
            else:
                if collect_stats:
                    fstats.failed += 1
                keep = False
                failing_filter = filter_name
                break

        if keep:
            cleaned.append(record)
        else:
            logger.debug(
                "Reaction %s dropped by filter %s",
                record.reaction_id or f"idx-{idx}",
                failing_filter,
            )

    stats.n_output = len(cleaned)
    logger.info("Cleaning finished: kept %d reactions", stats.n_output)
    return cleaned, stats


def clean_reactions(
    rxn_smiles_list: Iterable[ReactionRecord],
    filters: Optional[List[ReactionFilter]] = None,
    drop_failed_parse: bool = True,
    strict: bool = True,
) -> List[ReactionRecord]:
    """
    Parse and clean a list of reaction SMILES.

    This is the core entry point for the cleaning pipeline.

    Args:
        rxn_smiles_list:
            An iterable of ReactionRecord objects. If reactants/reagents/products
            are empty, they will be parsed from reaction_smiles.

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
    cleaned, _ = _clean_reactions_internal(
        rxn_smiles_list=rxn_smiles_list,
        filters=filters,
        drop_failed_parse=drop_failed_parse,
        strict=strict,
        collect_stats=False,
    )
    return cleaned


def clean_reactions_with_report(
    rxn_smiles_list: Iterable[ReactionRecord],
    filters: Optional[List[ReactionFilter]] = None,
    drop_failed_parse: bool = True,
    strict: bool = True,
) -> Tuple[List[ReactionRecord], CleaningStats]:
    """
    Run the cleaning pipeline and return both the cleaned reactions and stats.
    """
    return _clean_reactions_internal(
        rxn_smiles_list=rxn_smiles_list,
        filters=filters,
        drop_failed_parse=drop_failed_parse,
        strict=strict,
        collect_stats=True,
    )


def clean_and_canonicalize(
    rxn_smiles_list: Iterable[ReactionRecord],
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
            Iterable of ReactionRecords.

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

    logger.info("Canonicalizing %d reactions (isomeric=%s)", len(cleaned), isomeric)
    canon_records: List[ReactionRecord] = []
    for rec in cleaned:
        try:
            canon_records.append(canonicalize_reaction(rec, isomeric=isomeric))
        except Exception:
            logger.exception(
                "Canonicalization failed for reaction %s",
                rec.reaction_id or rec.reaction_smiles or "<unknown>",
            )
            raise

    return canon_records


def basic_cleaning_pipeline(
    rxn_smiles_list: Iterable[ReactionRecord],
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
