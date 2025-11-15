# chemrxn_cleaner/reporting.py

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Sequence, Tuple, Dict, Any

from .types import ReactionRecord


@dataclass
class CleaningReport:
    """
    Summary statistics of a reaction dataset before and after cleaning.
    """

    total_before: int
    total_after: int

    avg_n_reactants: float
    avg_n_reagents: float
    avg_n_products: float

    median_n_reactants: float
    median_n_reagents: float
    median_n_products: float

    def n_dropped(self) -> int:
        return self.total_before - self.total_after

    def drop_rate(self) -> float:
        if self.total_before == 0:
            return 0.0
        return self.n_dropped() / self.total_before

    def to_dict(self) -> dict:
        return {
            "total_before": self.total_before,
            "total_after": self.total_after,
            "n_dropped": self.n_dropped(),
            "drop_rate": self.drop_rate(),
            "avg_n_reactants": self.avg_n_reactants,
            "avg_n_reagents": self.avg_n_reagents,
            "avg_n_products": self.avg_n_products,
            "median_n_reactants": self.median_n_reactants,
            "median_n_reagents": self.median_n_reagents,
            "median_n_products": self.median_n_products,
        }

    def pretty_print(self) -> None:
        """
        Print a human-readable summary to stdout.
        """
        print("=" * 60)
        print(" ChemRxn-Cleaner: Cleaning Summary")
        print("=" * 60)
        print(f"  Total reactions before cleaning : {self.total_before}")
        print(f"  Total reactions after  cleaning : {self.total_after}")
        print(f"  Dropped reactions               : {self.n_dropped()}")
        print(f"  Drop rate                       : {self.drop_rate():.2%}")
        print("-" * 60)
        print("  Reactants count per reaction:")
        print(f"    avg    : {self.avg_n_reactants:.2f}")
        print(f"    median : {self.median_n_reactants:.2f}")
        print("  Reagents count per reaction:")
        print(f"    avg    : {self.avg_n_reagents:.2f}")
        print(f"    median : {self.median_n_reagents:.2f}")
        print("  Products count per reaction:")
        print(f"    avg    : {self.avg_n_products:.2f}")
        print(f"    median : {self.median_n_products:.2f}")
        print("=" * 60)


# ------------------------ helper functions ------------------------ #


def _median(values: List[int]) -> float:
    if not values:
        return 0.0
    values_sorted = sorted(values)
    n = len(values_sorted)
    mid = n // 2
    if n % 2 == 1:
        return float(values_sorted[mid])
    else:
        return (values_sorted[mid - 1] + values_sorted[mid]) / 2.0


# ------------------------ public API ------------------------ #


def summarize_cleaning(
    raw_reactions: Sequence[Tuple[str, Dict[str, Any]]],
    cleaned_reactions: Sequence[ReactionRecord],
) -> CleaningReport:
    """
    Compare raw vs cleaned reactions and produce a summary report.

    Note:
        This function currently only computes aggregate statistics.
        It does not break down drop reasons by individual filters.
        (That can be added later if we track per-filter decisions.)
    """
    total_before = len(raw_reactions)
    total_after = len(cleaned_reactions)

    if total_after == 0:
        # Avoid division by zero; return a mostly-empty report.
        return CleaningReport(
            total_before=total_before,
            total_after=total_after,
            avg_n_reactants=0.0,
            avg_n_reagents=0.0,
            avg_n_products=0.0,
            median_n_reactants=0.0,
            median_n_reagents=0.0,
            median_n_products=0.0,
        )

    n_reactants_list: List[int] = []
    n_reagents_list: List[int] = []
    n_products_list: List[int] = []

    for rec in cleaned_reactions:
        n_reactants_list.append(len(rec.reactants))
        n_reagents_list.append(len(rec.reagents))
        n_products_list.append(len(rec.products))

    avg_n_reactants = sum(n_reactants_list) / total_after
    avg_n_reagents = sum(n_reagents_list) / total_after
    avg_n_products = sum(n_products_list) / total_after

    median_n_reactants = _median(n_reactants_list)
    median_n_reagents = _median(n_reagents_list)
    median_n_products = _median(n_products_list)

    return CleaningReport(
        total_before=total_before,
        total_after=total_after,
        avg_n_reactants=avg_n_reactants,
        avg_n_reagents=avg_n_reagents,
        avg_n_products=avg_n_products,
        median_n_reactants=median_n_reactants,
        median_n_reagents=median_n_reagents,
        median_n_products=median_n_products,
    )
