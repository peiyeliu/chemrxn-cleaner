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


def _write_json(
    records: Sequence[ReactionRecord], out_path: Path, *, indent: int
) -> None:
    payload = [record.to_dict() for record in records]
    logger.info("Exporting %d reactions to JSON at %s", len(records), out_path)
    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=indent)


def _write_csv(records: Sequence[ReactionRecord], out_path: Path) -> None:
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


def export_reaction_records(
    records: Sequence[ReactionRecord],
    path: str | Path,
    *,
    fmt: str | None = None,
    indent: int = 2,
) -> None:
    """Serialize reactions to the provided path as JSON or CSV.

    Args:
        records: Reaction records to serialize.
        path: Destination path for the file.
        fmt: Output format; ``"json"`` or ``"csv"``. When omitted, the format is
            inferred from the file extension.
        indent: Indentation level for JSON output (ignored for CSV).
    """
    out_path = _ensure_path(path)
    resolved_fmt = (fmt or out_path.suffix.lstrip(".")).lower()
    if not resolved_fmt:
        raise ValueError("File format not provided; supply fmt or a .json/.csv path.")

    if resolved_fmt == "json":
        _write_json(records, out_path, indent=indent)
    elif resolved_fmt == "csv":
        _write_csv(records, out_path)
    else:
        raise ValueError(
            f"Unsupported export format '{resolved_fmt}'. Use 'json' or 'csv'."
        )
