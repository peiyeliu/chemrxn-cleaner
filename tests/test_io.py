import pytest

from chemrxn_cleaner.types import ReactionRecord
from chemrxn_cleaner.io import (
    export_reaction_records_to_json,
    load_reaction_records_from_json,
    export_reaction_records_to_csv,
    load_reaction_records_from_csv,
)


def _sample_records():
    return [
        ReactionRecord(
            raw="CCC.CCO>>CCOCC",
            reactants=["CCC", "CCO"],
            reagents=[],
            products=["CCOCC"],
            meta={"source": "unit-test", "index": 1},
        ),
        ReactionRecord(
            raw="O=O>>O",
            reactants=["O=O"],
            reagents=["[H][H]"],
            products=["O"],
            meta=None,
        ),
    ]


def test_json_round_trip(tmp_path):
    records = _sample_records()
    path = tmp_path / "records.json"

    export_reaction_records_to_json(records, path)
    loaded = load_reaction_records_from_json(path)

    assert loaded == records


def test_csv_round_trip(tmp_path):
    records = _sample_records()
    path = tmp_path / "records.csv"

    export_reaction_records_to_csv(records, path)
    loaded = load_reaction_records_from_csv(path)

    assert loaded == records


def test_csv_missing_columns(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("raw,reactants\nfoo,[]\n", encoding="utf-8")

    with pytest.raises(ValueError):
        load_reaction_records_from_csv(path)


def test_csv_invalid_json_list(tmp_path):
    path = tmp_path / "bad_list.csv"
    path.write_text(
        "raw,reactants,reagents,products,meta\n" 'foo,"not json",[],[],{}\n',
        encoding="utf-8",
    )

    with pytest.raises(ValueError):
        load_reaction_records_from_csv(path)
