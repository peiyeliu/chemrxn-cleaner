# chemrxn_cleaner/ml/utils.py
from __future__ import annotations

import random
from dataclasses import asdict
from typing import Iterable, List, Tuple

import pandas as pd

from chemrxn_cleaner.types import ReactionRecord


def records_to_dataframe(records: Iterable[ReactionRecord]) -> pd.DataFrame:
    """
    Convert a list of ReactionRecord objects into a pandas.DataFrame.

    - Keeps all fields (including extra_* dicts) by flattening dataclasses.
    - Good for quick EDA / exporting cleaned data.
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
    """
    Simple random split. You can later replace this with
    scaffold- / time-based splitting but keep the same interface.
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
        return [records[i] for i in idxs_]

    return pick(train_idx), pick(valid_idx), pick(test_idx)
