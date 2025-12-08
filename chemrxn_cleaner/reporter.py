# chemrxn_cleaner/reporter.py

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable


@dataclass
class FilterStats:
    """
    Stats for a specific filter
    """

    name: str
    applied: int = 0
    passed: int = 0
    failed: int = 0

    def format(self) -> str:
        """
        Human-friendly one-line summary.
        """
        return (
            f"{self.name}: applied={self.applied}, "
            f"passed={self.passed}, failed={self.failed}"
        )

    def __str__(self) -> str:
        return self.format()


@dataclass
class CleaningStats:
    """
    Stats for the overall Cleaning process
    """

    n_input: int = 0
    n_output: int = 0
    n_failed_parse: int = 0
    per_filter: Dict[str, FilterStats] = field(default_factory=dict)

    def _merge_filter_stats(self, source: Dict[str, FilterStats]) -> None:
        """
        Accumulate filter-level counters from another mapping.
        """
        for name, fstats in source.items():
            merged = self.per_filter.setdefault(name, FilterStats(name=name))
            merged.applied += fstats.applied
            merged.passed += fstats.passed
            merged.failed += fstats.failed

    def __iadd__(self, other: "CleaningStats") -> "CleaningStats":
        if not isinstance(other, CleaningStats):
            return NotImplemented

        self.n_input += other.n_input
        self.n_output += other.n_output
        self.n_failed_parse += other.n_failed_parse
        self._merge_filter_stats(other.per_filter)
        return self

    def __add__(self, other: "CleaningStats") -> "CleaningStats":
        if not isinstance(other, CleaningStats):
            return NotImplemented

        combined = CleaningStats(
            n_input=self.n_input,
            n_output=self.n_output,
            n_failed_parse=self.n_failed_parse,
        )
        combined._merge_filter_stats(self.per_filter)
        combined += other
        return combined

    @classmethod
    def combine(cls, stats_list: Iterable["CleaningStats"]) -> "CleaningStats":
        """
        Combine multiple CleaningStats objects (e.g., from parallel runs).
        """
        combined = cls()
        for stats in stats_list:
            combined += stats
        return combined

    def summary(self, include_filters: bool = True) -> str:
        """
        Return a formatted string summarizing counters.
        """
        lines = [
            f"Input: {self.n_input}, output: {self.n_output}, "
            f"failed_parse: {self.n_failed_parse}"
        ]
        if include_filters and self.per_filter:
            lines.append("Per-filter:")
            for name in sorted(self.per_filter.keys()):
                lines.append(f"  {self.per_filter[name]}")
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.summary()
