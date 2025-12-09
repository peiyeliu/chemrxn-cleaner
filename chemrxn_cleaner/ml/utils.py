# chemrxn_cleaner/ml/utils.py
from __future__ import annotations

import random
from dataclasses import asdict
from typing import Iterable, List, Tuple

import pandas as pd

from chemrxn_cleaner.types import ReactionRecord


def records_to_dataframe(records: Iterable[ReactionRecord]) -> pd.DataFrame:
    """Convert reaction records into a ``pandas.DataFrame``.

    Args:
        records: Iterable of ``ReactionRecord`` instances.

    Returns:
        DataFrame containing flattened record fields; list fields are joined
        with ``" | "`` separators for readability.
    """
    rows = []
    for r in records:
        row = asdict(r)
        # Flatten list fields to simple strings for CSV (optional)
        for key, value in list(row.items()):
            if isinstance(value, list):
                row[key] = " | ".join(value)
        rows.append(row)
    return pd.DataFrame(rows)


def train_valid_test_split(
    records: List[ReactionRecord],
    train_ratio: float = 0.8,
    valid_ratio: float = 0.1,
    seed: int = 0,
) -> Tuple[List[ReactionRecord], List[ReactionRecord], List[ReactionRecord]]:
    """Randomly split records into train/valid/test partitions.

    Args:
        records: Dataset to split.
        train_ratio: Fraction of examples allocated to the training set.
        valid_ratio: Fraction allocated to the validation set.
        seed: Random seed for deterministic shuffling.

    Returns:
        Tuple of ``(train, valid, test)`` record lists.
    """
    rng = random.Random(seed)
    idxs = list(range(len(records)))
    rng.shuffle(idxs)

    n = len(records)
    n_train = int(n * train_ratio)
    n_valid = int(n * valid_ratio)

    train_idx = idxs[:n_train]
    valid_idx = idxs[n_train : n_train + n_valid]
    test_idx = idxs[n_train + n_valid :]

    def pick(idxs_):
        """Select records by shuffled indices.

        Args:
            idxs_: Indices to extract from the records list.

        Returns:
            Ordered list of ``ReactionRecord`` instances corresponding to
            ``idxs_``.
        """
        return [records[i] for i in idxs_]

    return pick(train_idx), pick(valid_idx), pick(test_idx)
