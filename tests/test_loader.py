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


def test_extract_ord_reaction_smiles_procedure_yield(monkeypatch):
    class DummyProduct:
        def __init__(self, smiles: str, yield_value: float | None):
            self.smiles = smiles
            self.yield_value = yield_value

    class DummyOutcome:
        def __init__(self, products):
            self.products = products

    class DummyReaction:
        def __init__(self, reaction_id, smiles, outcomes, flat):
            self.reaction_id = reaction_id
            self.smiles = smiles
            self.outcomes = outcomes
            self.flat = flat

    dummy_dataset = types.SimpleNamespace(
        reactions=[
            DummyReaction(
                "rxn-42",
                "SMI",
                [
                    DummyOutcome(
                        [
                            DummyProduct("P1", 90.0),
                            DummyProduct("P2", None),
                        ]
                    )
                ],
                {
                    "setup.temperature": 298,
                    "conditions.time": 5,
                    "other.field": "skip",
                },
            )
        ]
    )

    def fake_load_message(path, message_cls):
        return dummy_dataset

    def fake_get_reaction_smiles(message, **kwargs):
        return message.smiles

    def fake_smiles_from_compound(product, canonical=True):
        return product.smiles

    def fake_get_product_yield(product, as_measurement=False):
        return product.yield_value

    def fake_message_to_row(rxn):
        return rxn.flat

    monkeypatch.setattr(loader, "load_message", fake_load_message)
    monkeypatch.setattr(loader, "get_reaction_smiles", fake_get_reaction_smiles)
    monkeypatch.setattr(loader, "smiles_from_compound", fake_smiles_from_compound)
    monkeypatch.setattr(loader, "get_product_yield", fake_get_product_yield)
    monkeypatch.setattr(loader, "message_to_row", fake_message_to_row)

    reactions = loader.extract_ord_reaction_smiles_procedure_yield("dummy.pb")

    assert reactions == [
        (
            "SMI",
            {
                "reaction_id": "rxn-42",
                "procedure": {
                    "setup.temperature": 298,
                    "conditions.time": 5,
                },
                "yields": [
                    {"product_smiles": "P1", "yield_percent": 90.0},
                ],
            },
        )
    ]
