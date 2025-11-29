import types
import json

import pytest

from chemrxn_cleaner.io import loader


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
    class DummyReaction:
        def __init__(self, reaction_id: str, smiles: str):
            self.reaction_id = reaction_id
            self.smiles = smiles

    dummy_dataset = types.SimpleNamespace(
        reactions=[
            DummyReaction("rxn-1", "A>B>C"),
            DummyReaction("", "D>E>F"),
        ]
    )

    def fake_load_message(path, message_cls):
        assert path == "dummy.pb"
        return dummy_dataset

    def fake_get_reaction_smiles(message, **kwargs):
        return message.smiles

    monkeypatch.setattr(loader, "load_message", fake_load_message)
    monkeypatch.setattr(loader, "get_reaction_smiles", fake_get_reaction_smiles)

    reactions = loader.load_ord("dummy.pb")

    assert reactions == [
        ("A>B>C", {"reaction_id": "rxn-1", "reaction_index": 0}),
        ("D>E>F", {"reaction_index": 1}),
    ]


def test_load_ord_applies_meta_extractor(monkeypatch):
    class DummyReaction:
        def __init__(self, reaction_id: str, smiles: str):
            self.reaction_id = reaction_id
            self.smiles = smiles

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

    monkeypatch.setattr(loader, "load_message", fake_load_message)
    monkeypatch.setattr(loader, "get_reaction_smiles", fake_get_reaction_smiles)

    def meta_extractor(reaction):
        return {"extra": f"meta-{reaction.smiles}"}

    reactions = loader.load_ord("dummy.pb", meta_extractor=meta_extractor)

    assert reactions == [
        ("A>B>C", {"reaction_id": "rxn-1", "reaction_index": 0, "extra": "meta-A>B>C"}),
        ("D>E>F", {"reaction_index": 1, "extra": "meta-D>E>F"}),
    ]
