import types

import pytest

from chemrxn_cleaner import loader


def test_load_uspto_rsmi_without_meta(tmp_path):
    data = "CCO.CC>O>CO\nCCO>O>CC\tfoo\tbar\n"
    path = tmp_path / "sample.rsmi"
    path.write_text(data, encoding="utf-8")

    reactions = loader.load_uspto_rsmi(str(path), keep_meta=False)

    assert reactions == [
        ("CCO.CC>O>CO", {}),
        ("CCO>O>CC", {}),
    ]


def test_load_uspto_rsmi_with_meta(tmp_path):
    data = "CCO>O>CO\tfield1\tfield2\n"
    path = tmp_path / "sample_meta.rsmi"
    path.write_text(data, encoding="utf-8")

    reactions = loader.load_uspto_rsmi(str(path), keep_meta=True)

    assert reactions == [
        ("CCO>O>CO", {"fields": ["field1", "field2"]}),
    ]


def test_load_ord_pb_reaction_smiles_returns_meta(monkeypatch):
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

    reactions = loader.load_ord_pb_reaction_smiles("dummy.pb")

    assert reactions == [
        ("A>B>C", {"reaction_id": "rxn-1", "reaction_index": 0}),
        ("D>E>F", {"reaction_index": 1}),
    ]


def test_load_ord_pb_reaction_smiles_applies_meta_extractor(monkeypatch):
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

    reactions = loader.load_ord_pb_reaction_smiles(
        "dummy.pb", meta_extractor=meta_extractor
    )

    assert reactions == [
        ("A>B>C", {"reaction_id": "rxn-1", "reaction_index": 0, "extra": "meta-A>B>C"}),
        ("D>E>F", {"reaction_index": 1, "extra": "meta-D>E>F"}),
    ]
