import types
import json

import pytest

from chemrxn_cleaner.io import loader
from chemrxn_cleaner.types import ReactionRecord, YieldType


def test_load_uspto_without_meta(tmp_path):
    data = "CCO.CC>O>CO\nCCO>O>CC\tfoo\tbar\n"
    path = tmp_path / "sample.rsmi"
    path.write_text(data, encoding="utf-8")

    reactions = loader.load_uspto(str(path), keep_meta=False)

    assert reactions == [
        ("CCO.CC>O>CO", {}),
        ("CCO>O>CC", {}),
    ]


def test_load_uspto_with_meta(tmp_path):
    data = "CCO>O>CO\tfield1\tfield2\n"
    path = tmp_path / "sample_meta.rsmi"
    path.write_text(data, encoding="utf-8")

    reactions = loader.load_uspto(str(path), keep_meta=True)

    assert reactions == [
        ("CCO>O>CO", {"fields": ["field1", "field2"]}),
    ]


def test_load_csv_infers_meta(tmp_path):
    data = (
        "reactant1,reactant2,reagent,product,temp,note\n"
        "CCO,CCBr,NaOH,CCOCC,25C,batch1\n"
        "CCN,,,CCNH2,,batch2\n"
    )
    path = tmp_path / "rxns.csv"
    path.write_text(data, encoding="utf-8")

    reactions = loader.load_csv(
        str(path),
        reactant_columns=["reactant1", "reactant2"],
        reagent_columns=["reagent"],
        product_columns=["product"],
    )

    assert reactions == [
        ("CCO.CCBr>NaOH>CCOCC", {"temp": "25C", "note": "batch1"}),
        ("CCN>>CCNH2", {"temp": "", "note": "batch2"}),
    ]


def test_load_csv_custom_meta(tmp_path):
    data = "reactant,reagent,product,ref,source\n" "A,,C,foo,bar\n"
    path = tmp_path / "rxn_single.csv"
    path.write_text(data, encoding="utf-8")

    reactions = loader.load_csv(
        str(path),
        reactant_columns=["reactant"],
        product_columns=["product"],
        meta_columns=["ref"],
    )

    assert reactions == [("A>>C", {"ref": "foo"})]


def test_load_csv_combined_column(tmp_path):
    data = "rxn_smiles,tag\n" "A.B>C>D,batch\n" ",empty\n"
    path = tmp_path / "rxn_combined.csv"
    path.write_text(data, encoding="utf-8")

    reactions = loader.load_csv(
        str(path),
        reaction_smiles_column="rxn_smiles",
    )

    assert reactions == [("A.B>C>D", {"tag": "batch"})]


def test_load_json_with_mapper(tmp_path):
    payload = [
        {"reactants": "A.B", "products": "C", "meta": {"id": 1}},
        {"skip": True},
    ]
    path = tmp_path / "rxns.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    def mapper(entry):
        if entry.get("skip"):
            return None
        return (
            f"{entry['reactants']}>>{entry['products']}",
            {"id": entry["meta"]["id"]},
        )

    reactions = loader.load_json(str(path), mapper)

    assert reactions == [("A.B>>C", {"id": 1})]


def test_load_ord_returns_meta(monkeypatch):
    class DummyProduct:
        def __init__(self, yield_value: float | None):
            self.yield_value = yield_value

    class DummyOutcome:
        def __init__(self, products):
            self.products = products

    class DummyReaction:
        def __init__(self, reaction_id: str, smiles: str, outcomes, flat):
            self.reaction_id = reaction_id
            self.smiles = smiles
            self.outcomes = outcomes
            self.flat = flat

    dummy_dataset = types.SimpleNamespace(
        reactions=[
            DummyReaction(
                "rxn-1",
                "A>B>C",
                [DummyOutcome([DummyProduct(82.5)])],
                {"conditions.temperature": 25, "conditions.time": 2},
            ),
            DummyReaction(
                "",
                "D>E>F",
                [DummyOutcome([DummyProduct(None)])],
                {},
            ),
        ]
    )

    def fake_load_message(path, message_cls):
        assert path == "dummy.pb"
        return dummy_dataset

    def fake_get_reaction_smiles(message, **kwargs):
        return message.smiles

    def fake_get_product_yield(product, as_measurement=False):
        return product.yield_value

    def fake_message_to_row(reaction):
        return reaction.flat

    monkeypatch.setattr(loader, "load_message", fake_load_message)
    monkeypatch.setattr(loader, "get_reaction_smiles", fake_get_reaction_smiles)
    monkeypatch.setattr(loader, "get_product_yield", fake_get_product_yield)
    monkeypatch.setattr(loader, "message_to_row", fake_message_to_row)

    reactions = loader.load_ord("dummy.pb")

    assert isinstance(reactions[0], ReactionRecord)
    assert reactions[0].reaction_id == "rxn-1"
    assert reactions[0].reactants == ["A"]
    assert reactions[0].products == ["C"]
    assert reactions[0].yield_value == 82.5
    assert reactions[0].yield_type == YieldType.OTHER
    assert reactions[0].temperature_c == 25
    assert reactions[0].time_hours == 2
    assert reactions[0].extra_metadata["reaction_index"] == 0

    assert reactions[1].reaction_id == ""
    assert reactions[1].yield_value is None
    assert reactions[1].yield_type == YieldType.NONE
    assert reactions[1].extra_metadata["reaction_index"] == 1


def test_load_ord_applies_meta_extractor(monkeypatch):
    class DummyReaction:
        def __init__(self, reaction_id: str, smiles: str):
            self.reaction_id = reaction_id
            self.smiles = smiles
            self.outcomes = []

    dummy_dataset = types.SimpleNamespace(
        reactions=[
            DummyReaction("rxn-1", "A>B>C"),
            DummyReaction("", "D>E>F"),
        ]
    )

    def fake_load_message(path, message_cls):
        return dummy_dataset

    def fake_get_reaction_smiles(message, **kwargs):
        return message.smiles

    def fake_get_product_yield(product, as_measurement=False):
        return None

    def fake_message_to_row(reaction):
        return {}

    monkeypatch.setattr(loader, "load_message", fake_load_message)
    monkeypatch.setattr(loader, "get_reaction_smiles", fake_get_reaction_smiles)
    monkeypatch.setattr(loader, "get_product_yield", fake_get_product_yield)
    monkeypatch.setattr(loader, "message_to_row", fake_message_to_row)

    def meta_extractor(reaction):
        return {"extra": f"meta-{reaction.smiles}"}

    reactions = loader.load_ord("dummy.pb", meta_extractor=meta_extractor)

    assert reactions[0].extra_metadata["extra"] == "meta-A>B>C"
    assert reactions[0].reaction_id == "rxn-1"
    assert reactions[1].extra_metadata["extra"] == "meta-D>E>F"
    assert reactions[1].reaction_id == ""
