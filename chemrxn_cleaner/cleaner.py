# chemrxn_cleaner/cleaner.py

from __future__ import annotations

import logging
from typing import Iterable, List, Optional, Tuple

from .filters import ReactionFilter
from .parser import canonicalize_reaction, parse_reaction_smiles
from .reporter import CleaningStats, FilterStats
from .types import ReactionRecord

logger = logging.getLogger(__name__)


def _get_filter_name(f: ReactionFilter) -> str:
    """Return a human-readable name for a reaction filter.

    Args:
        f: Filter callable used in the cleaning pipeline.

    Returns:
        The callable ``__name__`` when present, otherwise the class name.
    """
    return getattr(f, "__name__", f.__class__.__name__)


def _ensure_filter_stats(stats: CleaningStats, filter_name: str) -> FilterStats:
    """Ensure the per-filter stats entry exists and return it.

    Args:
        stats: Aggregate cleaning stats to populate.
        filter_name: Name of the filter being evaluated.

    Returns:
        The ``FilterStats`` instance corresponding to ``filter_name``.
    """
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
    """Core implementation for the cleaning pipeline.

    Args:
        rxn_smiles_list: Incoming iterable of reaction records to clean.
        filters: Optional list of ReactionFilter callables to apply. When
            ``None``, no filters are executed; pass ``default_filters()`` to
            apply the recommended stack.
        drop_failed_parse: Whether to drop reactions that fail SMILES parsing
            instead of propagating the error.
        strict: When True, enforce three '>' parts during SMILES parsing.
        collect_stats: When True, populate and return cleaning statistics.

    Returns:
        A tuple of ``(cleaned_records, stats)`` where ``stats`` always reflects
        the processing that occurred, even when ``collect_stats`` is False.

    Raises:
        Exception: Propagates parsing or canonicalization errors when
            ``drop_failed_parse`` is False.
    """
    filters = [] if filters is None else list(filters)

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
    """Parse, validate, and filter reaction records.

    Args:
        rxn_smiles_list: Iterable of reactions. Empty reactant/reagent/product
            fields are parsed from ``reaction_smiles`` when present.
        filters: Optional list of predicate callables. When omitted, no filters
            are applied; pass ``default_filters()`` for the recommended stack.
        drop_failed_parse: When True, quietly drop reactions that cannot be
            parsed; when False, propagate the parsing error.
        strict: Passed to ``parse_reaction_smiles``; enforces three '>' parts
            when True, otherwise pads missing parts.

    Returns:
        Cleaned reactions that passed all filters.
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
    """Run the cleaning pipeline and return results plus statistics.

    Args:
        rxn_smiles_list: Iterable of input reactions to process.
        filters: Optional list of filters to apply. When omitted, no filters
            are applied; use ``default_filters()`` for the recommended set.
        drop_failed_parse: Whether to drop reactions that cannot be parsed
            instead of raising.
        strict: When True, require reaction SMILES to contain three sections.

    Returns:
        A tuple of ``(cleaned_reactions, cleaning_stats)``.
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
    """Clean reactions and canonicalize all SMILES in one pass.

    Args:
        rxn_smiles_list: Iterable of ReactionRecords to clean.
        filters: Optional list of filter predicates. When None, no filters are
            run; pass ``default_filters()`` to apply the recommended stack.
        drop_failed_parse: Whether to drop reactions that fail parsing rather
            than raising.
        strict: Enforce exactly three ``>`` parts in reaction SMILES when True.
        isomeric: Preserve isomeric information during canonicalization.

    Returns:
        Reaction records with canonicalized reactants, reagents, and products.
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
