from __future__ import annotations

from pathlib import Path
from typing import List, Optional

from chemrxn_cleaner.io import loader
from chemrxn_cleaner.parser import parse_reaction_smiles
from chemrxn_cleaner.types import ReactionRecord, YieldType

RESOURCES_DIR = Path(__file__).parent / "resources"


def _assert_reaction_records(reactions: List[ReactionRecord], label: str) -> None:
    """Confirm that a loader returned at least one ReactionRecord."""
    assert reactions, f"{label} loader returned no reactions."
    assert all(isinstance(rec, ReactionRecord) for rec in reactions)


def _parse_yield(value: Optional[str]) -> Optional[float]:
    """Parse a yield value that may include a trailing '%'."""
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    text = text.rstrip("%")
    try:
        return float(text)
    except ValueError:
        return None


def test_load_uspto_sample() -> None:
    uspto_path = RESOURCES_DIR / "uspto_dataset-sample-small.rsmi"
    reactions = loader.load_uspto(str(uspto_path), keep_meta=True)
    _assert_reaction_records(reactions, "USPTO")
    assert reactions[0].source == "uspto"
    assert ">" in reactions[0].reaction_smiles
    assert reactions[0].extra_metadata["fields"][0] == "US03930836"


def test_load_csv_sample() -> None:
    csv_path = RESOURCES_DIR / "csv_dataset-sample.csv"

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

    reactions = loader.load_csv(
        str(csv_path),
        reaction_smiles_column="CanonicalizedReaction",
        delimiter="\t",
        skip_lines=2,  # Skip the comment lines at the top of the sample file.
        mapper=mapper,
    )
    _assert_reaction_records(reactions, "CSV")
    assert reactions[0].source == "csv"
    assert ">" in reactions[0].reaction_smiles


def test_load_json_sample() -> None:
    json_path = RESOURCES_DIR / "json_dataset-sample.json"

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

    reactions = loader.load_json(str(json_path), mapper)
    _assert_reaction_records(reactions, "JSON")
    assert reactions[0].source == "json"
    assert ">" in reactions[0].reaction_smiles


def test_load_ord_sample() -> None:
    ord_path = RESOURCES_DIR / "ord_dataset-sample.pb.gz"
    reactions = loader.load_ord(str(ord_path))
    _assert_reaction_records(reactions, "ORD")
    assert reactions[0].source == "ord"
    assert ">" in reactions[0].reaction_smiles
