# chemrxn_cleaner/io.py
"""Utilities for serializing :class:`ReactionRecord` objects."""

from __future__ import annotations

from pathlib import Path
import json
import csv
from typing import Iterable, List, Sequence

from .types import ReactionRecord


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


def load_reaction_records_from_json(path: str | Path) -> List[ReactionRecord]:
    """Load :class:`ReactionRecord` objects from a JSON file."""
    in_path = _ensure_path(path)
    with in_path.open("r", encoding="utf-8") as fh:
        payload = json.load(fh)

    if not isinstance(payload, list):
        raise ValueError("JSON payload must be a list of reaction records")

    records: List[ReactionRecord] = []
    for idx, entry in enumerate(payload):
        if not isinstance(entry, dict):
            raise ValueError(
                f"Reaction entry at index {idx} must be a mapping, got {type(entry)!r}"
            )
        records.append(ReactionRecord.from_dict(entry))

    return records


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


def load_reaction_records_from_csv(path: str | Path) -> List[ReactionRecord]:
    """Load :class:`ReactionRecord` objects from a CSV file."""
    in_path = _ensure_path(path)
    with in_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        missing = {"raw", "reactants", "reagents", "products", "meta"} - set(
            reader.fieldnames or []
        )
        if missing:
            missing_list = ", ".join(sorted(missing))
            raise ValueError(f"CSV is missing required columns: {missing_list}")

        records: List[ReactionRecord] = []
        for row_idx, row in enumerate(reader):
            reactants = _parse_json_list(row.get("reactants"), row_idx, "reactants")
            reagents = _parse_json_list(row.get("reagents"), row_idx, "reagents")
            products = _parse_json_list(row.get("products"), row_idx, "products")
            meta = _parse_json_meta(row.get("meta"), row_idx)
            records.append(
                ReactionRecord(
                    raw=row.get("raw", "") or "",
                    reactants=reactants,
                    reagents=reagents,
                    products=products,
                    meta=meta,
                )
            )

    return records


def _parse_json_list(value: str | None, row_idx: int, field: str) -> List[str]:
    if not value:
        return []
    try:
        payload = json.loads(value)
    except json.JSONDecodeError as exc:  # pragma: no cover - error path
        raise ValueError(
            f"Row {row_idx}: unable to parse {field} as JSON list"
        ) from exc
    if not isinstance(payload, list):
        raise ValueError(
            f"Row {row_idx}: expected {field} to decode to a list, got {type(payload)!r}"
        )
    return [str(item) for item in payload]


def _parse_json_meta(value: str | None, row_idx: int) -> dict | None:
    if not value:
        return None
    try:
        payload = json.loads(value)
    except json.JSONDecodeError as exc:  # pragma: no cover - error path
        raise ValueError(f"Row {row_idx}: unable to parse meta JSON") from exc
    if not isinstance(payload, dict):
        raise ValueError(
            f"Row {row_idx}: expected meta to decode to a dict, got {type(payload)!r}"
        )
    return payload
