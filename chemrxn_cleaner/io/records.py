# chemrxn_cleaner/io/records.py
"""Utilities for serializing :class:`ReactionRecord` objects."""

from __future__ import annotations

from pathlib import Path
import json
import csv
from typing import Iterable, Sequence

from ..types import ReactionRecord


def _ensure_path(path: str | Path) -> Path:
    return path if isinstance(path, Path) else Path(path)


def export_reaction_records_to_json(
    records: Sequence[ReactionRecord],
    path: str | Path,
    *,
    indent: int = 2,
) -> None:
    """Write the provided reactions to ``path`` as a JSON list."""
    out_path = _ensure_path(path)
    payload = [record.to_dict() for record in records]
    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=indent)


def export_reaction_records_to_csv(
    records: Sequence[ReactionRecord],
    path: str | Path,
) -> None:
    """Write reactions to ``path`` in CSV format."""
    out_path = _ensure_path(path)
    fieldnames = ["raw", "reactants", "reagents", "products", "meta"]

    def dumps(obj: Iterable[str] | dict | None) -> str:
        if obj is None:
            return ""
        return json.dumps(obj, ensure_ascii=False)

    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "raw": record.raw,
                    "reactants": dumps(list(record.reactants)),
                    "reagents": dumps(list(record.reagents)),
                    "products": dumps(list(record.products)),
                    "meta": dumps(record.meta),
                }
            )
