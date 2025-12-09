# chemrxn_cleaner/io/writer.py
"""Utilities for serializing :class:`ReactionRecord` objects."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Iterable, Sequence

from ..types import ReactionRecord

logger = logging.getLogger(__name__)


def _ensure_path(path: str | Path) -> Path:
    """Convert a string or Path-like object into a ``Path`` instance.

    Args:
        path: Target filesystem path.

    Returns:
        Normalized ``Path`` object.
    """
    return path if isinstance(path, Path) else Path(path)


def export_reaction_records_to_json(
    records: Sequence[ReactionRecord],
    path: str | Path,
    *,
    indent: int = 2,
) -> None:
    """Write the provided reactions to ``path`` as a JSON list.

    Args:
        records: Reaction records to serialize.
        path: Destination path for the JSON file.
        indent: Indentation level for ``json.dump``.
    """
    out_path = _ensure_path(path)
    payload = [record.to_dict() for record in records]
    logger.info("Exporting %d reactions to JSON at %s", len(records), out_path)
    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=indent)


def export_reaction_records_to_csv(
    records: Sequence[ReactionRecord],
    path: str | Path,
) -> None:
    """Write reactions to ``path`` in CSV format.

    Args:
        records: Reaction records to export.
        path: Destination path for the CSV file.
    """
    out_path = _ensure_path(path)
    fieldnames = [
        "reaction_smiles",
        "reactants",
        "reagents",
        "products",
        "extra_metadata",
    ]

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
                    "reaction_smiles": record.reaction_smiles,
                    "reactants": dumps(list(record.reactants)),
                    "reagents": dumps(list(record.reagents)),
                    "products": dumps(list(record.products)),
                    "extra_metadata": dumps(record.extra_metadata),
                }
            )
    logger.info("Exported %d reactions to CSV at %s", len(records), out_path)
