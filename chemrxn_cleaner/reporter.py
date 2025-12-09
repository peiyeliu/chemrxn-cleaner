# chemrxn_cleaner/reporter.py

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable


@dataclass
class FilterStats:
    """Accumulates counts for an individual reaction filter.

    Attributes:
        name: Human-readable filter identifier.
        applied: Number of times the filter was evaluated.
        passed: Number of reactions accepted by the filter.
        failed: Number of reactions rejected by the filter.
    """

    name: str
    applied: int = 0
    passed: int = 0
    failed: int = 0

    def format(self) -> str:
        """Return a concise one-line summary of the counts."""
        return (
            f"{self.name}: applied={self.applied}, "
            f"passed={self.passed}, failed={self.failed}"
        )

    def __str__(self) -> str:
        return self.format()


@dataclass
class CleaningStats:
    """Aggregate statistics for the overall cleaning process.

    Attributes:
        n_input: Total number of reactions processed.
        n_output: Total number of reactions kept after filtering.
        n_failed_parse: Number of reactions that failed SMILES parsing.
        per_filter: Mapping of filter name to per-filter counters.
    """

    n_input: int = 0
    n_output: int = 0
    n_failed_parse: int = 0
    per_filter: Dict[str, FilterStats] = field(default_factory=dict)

    def _merge_filter_stats(self, source: Dict[str, FilterStats]) -> None:
        """Accumulate filter-level counters from another mapping.

        Args:
            source: Mapping of filter names to stats to merge into this object.
        """
        for name, fstats in source.items():
            merged = self.per_filter.setdefault(name, FilterStats(name=name))
            merged.applied += fstats.applied
            merged.passed += fstats.passed
            merged.failed += fstats.failed

    def __iadd__(self, other: "CleaningStats") -> "CleaningStats":
        """In-place addition of counters from another ``CleaningStats``.

        Args:
            other: Cleaning statistics to combine into this instance.

        Returns:
            This instance after accumulation.
        """
        if not isinstance(other, CleaningStats):
            return NotImplemented

        self.n_input += other.n_input
        self.n_output += other.n_output
        self.n_failed_parse += other.n_failed_parse
        self._merge_filter_stats(other.per_filter)
        return self

    def __add__(self, other: "CleaningStats") -> "CleaningStats":
        """Return a new ``CleaningStats`` combining two instances.

        Args:
            other: Cleaning statistics to add to this instance.

        Returns:
            A new ``CleaningStats`` containing summed counters.
        """
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
        """Combine multiple ``CleaningStats`` objects (e.g., from parallel runs).

        Args:
            stats_list: Iterable of statistics objects to aggregate.

        Returns:
            Aggregate statistics across all provided objects.
        """
        combined = cls()
        for stats in stats_list:
            combined += stats
        return combined

    def summary(self, include_filters: bool = True) -> str:
        """Return a formatted string summarizing counters.

        Args:
            include_filters: Whether to include per-filter breakdowns.

        Returns:
            Multiline string representation of the aggregated counters.
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
