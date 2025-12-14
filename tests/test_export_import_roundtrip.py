from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Optional

from chemrxn_cleaner.io import export_reaction_records, loader
from chemrxn_cleaner.parser import parse_reaction_smiles
from chemrxn_cleaner.types import ReactionRecord, YieldType

RESOURCES_DIR = Path(__file__).parent / "resources"


def _parse_yield(value: Optional[str]) -> Optional[float]:
    """Parse a yield value that may include a trailing '%'."""
    if value is None:
        return None
    text = str(value).strip().rstrip("%")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _load_uspto() -> list[ReactionRecord]:
    return loader.load_uspto(
        str(RESOURCES_DIR / "uspto_dataset-sample-small.rsmi"), keep_meta=True
    )


def _load_csv() -> list[ReactionRecord]:
    def mapper(record: ReactionRecord, row: dict) -> ReactionRecord:
        record.reaction_id = (row.get("PatentNumber") or "").strip()
        record.extra_metadata["year"] = row.get("Year")
        yield_value = _parse_yield(
            row.get("TextMinedYield") or row.get("CalculatedYield")
        )
        if yield_value is not None:
            record.yield_value = yield_value
            record.yield_type = YieldType.OTHER
        return record

    return loader.load_csv(
        str(RESOURCES_DIR / "csv_dataset-sample.csv"),
        reaction_smiles_column="CanonicalizedReaction",
        delimiter="\t",
        skip_lines=2,
        mapper=mapper,
    )


def _load_json() -> list[ReactionRecord]:
    def mapper(entry: dict) -> Optional[ReactionRecord]:
        product_smiles = (entry.get("product") or {}).get("smiles")
        reactant_smiles = [
            item.get("smiles")
            for item in entry.get("reactants", [])
            if item.get("smiles")
        ]
        if not product_smiles or not reactant_smiles:
            return None

        reagents: list[str] = []
        base_smiles = (entry.get("base") or {}).get("smiles")
        if base_smiles:
            reagents.append(base_smiles)

        rxn_smiles = (
            f"{'.'.join(reactant_smiles)}>{'.'.join(reagents)}>{product_smiles}"
        )
        record = parse_reaction_smiles(rxn_smiles, strict=False)
        record.reaction_id = str(entry.get("Id", ""))
        record.source = "json"
        if entry.get("yield") not in (None, ""):
            try:
                record.yield_value = float(entry["yield"])
                record.yield_type = YieldType.OTHER
            except (TypeError, ValueError):
                pass
        return record

    return loader.load_json(
        str(RESOURCES_DIR / "json_dataset-sample.json"),
        mapper,
    )


def _load_ord() -> list[ReactionRecord]:
    return loader.load_ord(str(RESOURCES_DIR / "ord_dataset-sample.pb.gz"))


def _structures_equal(a: object, b: object) -> bool:
    """Recursively compare two JSON-like structures treating NaN as equal."""
    if isinstance(a, float) and isinstance(b, float):
        if math.isnan(a) and math.isnan(b):
            return True
    if isinstance(a, dict) and isinstance(b, dict):
        if a.keys() != b.keys():
            return False
        return all(_structures_equal(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return False
        return all(_structures_equal(x, y) for x, y in zip(a, b))
    return a == b


def test_export_import_roundtrip(tmp_path) -> None:
    loaders = {
        "uspto": _load_uspto,
        "csv": _load_csv,
        "json": _load_json,
        "ord": _load_ord,
    }

    for name, load_fn in loaders.items():
        original_records = load_fn()
        assert original_records, f"{name} loader returned no records"

        export_path = tmp_path / f"{name}_export.json"
        export_reaction_records(original_records, export_path, fmt="json")

        loaded_payload = json.loads(export_path.read_text())
        roundtrip_records = [ReactionRecord.from_dict(d) for d in loaded_payload]

        assert len(roundtrip_records) == len(original_records)
        for before, after in zip(original_records, roundtrip_records):
            assert _structures_equal(
                before.to_dict(), after.to_dict()
            ), f"{name} record did not roundtrip cleanly"
