# chemrxn_cleaner/reporter.py

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict


@dataclass
class FilterStats:
    """
    Stats for a specific filter
    """

    name: str
    applied: int = 0
    passed: int = 0
    failed: int = 0


@dataclass
class CleaningStats:
    """
    Stats for the overall Cleaning process
    """

    n_input: int = 0
    n_output: int = 0
    n_failed_parse: int = 0
    per_filter: Dict[str, FilterStats] = field(default_factory=dict)
