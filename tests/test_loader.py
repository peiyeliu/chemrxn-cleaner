import json
import types

import pytest

from chemrxn_cleaner.io import loader
from chemrxn_cleaner.parser import parse_reaction_smiles
from chemrxn_cleaner.types import ReactionRecord, YieldType


def test_load_uspto_without_meta(tmp_path):
    data = "CCO.CC>O>CO\nCCO>O>CC\tfoo\tbar\n"
    path = tmp_path / "sample.rsmi"
    path.write_text(data, encoding="utf-8")

    reactions = loader.load_uspto(str(path), keep_meta=False)

    assert len(reactions) == 2
    assert isinstance(reactions[0], ReactionRecord)
    assert reactions[0].reaction_smiles == "CCO.CC>O>CO"
    assert reactions[0].reactants == ["CCO", "CC"]
    assert reactions[0].reagents == ["O"]
    assert reactions[0].products == ["CO"]
    assert reactions[0].source == "uspto"
    assert reactions[0].source_file_path == str(path)
    assert reactions[0].extra_metadata == {}
    assert reactions[1].reaction_smiles == "CCO>O>CC"


def test_load_uspto_with_meta(tmp_path):
    data = "CCO>O>CO\tfield1\tfield2\n"
    path = tmp_path / "sample_meta.rsmi"
    path.write_text(data, encoding="utf-8")

    reactions = loader.load_uspto(str(path), keep_meta=True)

    assert len(reactions) == 1
    rec = reactions[0]
    assert isinstance(rec, ReactionRecord)
    assert rec.reaction_smiles == "CCO>O>CO"
    assert rec.extra_metadata["fields"] == ["field1", "field2"]
    assert rec.source == "uspto"
    assert rec.source_file_path == str(path)


def test_load_csv_infers_meta(tmp_path):
    data = (
        "reactant1,reactant2,reagent,product,temp,note\n"
        "CCO,CCBr,NaOH,CCOCC,25C,batch1\n"
        "CCN,,,CCNH2,,batch2\n"
    )
    path = tmp_path / "rxns.csv"
    path.write_text(data, encoding="utf-8")

    def mapper(record, row):
        record.extra_metadata["temp"] = row["temp"]
        record.extra_metadata["note"] = row["note"]
        record.reaction_id = row["note"]
        return record

    reactions = loader.load_csv(
        str(path),
        reactant_columns=["reactant1", "reactant2"],
        reagent_columns=["reagent"],
        product_columns=["product"],
        mapper=mapper,
    )

    assert len(reactions) == 2
    assert isinstance(reactions[0], ReactionRecord)
    assert reactions[0].reaction_smiles == "CCO.CCBr>NaOH>CCOCC"
    assert reactions[0].reactants == ["CCO", "CCBr"]
    assert reactions[0].reagents == ["NaOH"]
    assert reactions[0].products == ["CCOCC"]
    assert reactions[0].reaction_id == "batch1"
    assert reactions[0].extra_metadata == {"temp": "25C", "note": "batch1"}
    assert reactions[1].reaction_smiles == "CCN>>CCNH2"
    assert reactions[1].reactants == ["CCN"]
    assert reactions[1].reagents == []
    assert reactions[1].products == ["CCNH2"]
    assert reactions[1].extra_metadata == {"temp": "", "note": "batch2"}


def test_load_csv_custom_meta(tmp_path):
    data = "reactant,reagent,product,ref,source\nA,,C,foo,bar\n"
    path = tmp_path / "rxn_single.csv"
    path.write_text(data, encoding="utf-8")

    def mapper(record, row):
        record.reaction_id = row["ref"]
        record.extra_metadata["source"] = row["source"]
        return record

    reactions = loader.load_csv(
        str(path),
        reactant_columns=["reactant"],
        product_columns=["product"],
        mapper=mapper,
    )

    assert len(reactions) == 1
    assert isinstance(reactions[0], ReactionRecord)
    assert reactions[0].reaction_smiles == "A>>C"
    assert reactions[0].reaction_id == "foo"
    assert reactions[0].extra_metadata["source"] == "bar"


def test_load_csv_combined_column(tmp_path):
    data = "rxn_smiles,tag\nA.B>C>D,batch\n,empty\n"
    path = tmp_path / "rxn_combined.csv"
    path.write_text(data, encoding="utf-8")

    def mapper(record, row):
        record.extra_metadata["tag"] = row["tag"]
        record.success = row["tag"] == "batch"
        return record

    reactions = loader.load_csv(
        str(path),
        reaction_smiles_column="rxn_smiles",
        mapper=mapper,
    )

    assert len(reactions) == 1
    assert reactions[0].reaction_smiles == "A.B>C>D"
    assert reactions[0].extra_metadata["tag"] == "batch"
    assert reactions[0].success is True


def test_load_csv_skips_initial_lines(tmp_path):
    data = "#comment header\n#skip me too\nreactant,product\nA,C\nB,D\n"
    path = tmp_path / "rxn_skip.csv"
    path.write_text(data, encoding="utf-8")

    reactions = loader.load_csv(
        str(path),
        reactant_columns=["reactant"],
        product_columns=["product"],
        skip_lines=2,
    )

    assert len(reactions) == 2
    assert reactions[0].reaction_smiles == "A>>C"
    assert reactions[1].reaction_smiles == "B>>D"


def test_load_csv_rejects_negative_skip(tmp_path):
    path = tmp_path / "rxn_skip_negative.csv"
    path.write_text("reactant,product\nA,C\n", encoding="utf-8")

    with pytest.raises(ValueError):
        loader.load_csv(
            str(path),
            reactant_columns=["reactant"],
            product_columns=["product"],
            skip_lines=-1,
        )


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
        record = parse_reaction_smiles(
            f"{entry['reactants']}>>{entry['products']}", strict=False
        )
        record.reaction_id = f"rxn-{entry['meta']['id']}"
        record.extra_metadata["id"] = entry["meta"]["id"]
        record.source = "json"
        return record

    reactions = loader.load_json(str(path), mapper)

    assert len(reactions) == 1
    assert isinstance(reactions[0], ReactionRecord)
    assert reactions[0].reaction_id == "rxn-1"
    assert reactions[0].reaction_smiles == "A.B>>C"
    assert reactions[0].reactants == ["A", "B"]
    assert reactions[0].products == ["C"]
    assert reactions[0].extra_metadata["id"] == 1
    assert reactions[0].source == "json"


def test_load_csv_strip_atom_mapping(tmp_path):
    data = "rxn_smiles\n[CH3:1].[OH-:2]>[Na+]>[CH3OH:1]\n"
    path = tmp_path / "rxn_map.csv"
    path.write_text(data, encoding="utf-8")

    reactions = loader.load_csv(
        str(path),
        reaction_smiles_column="rxn_smiles",
        strip_atom_mapping=True,
    )

    assert len(reactions) == 1
    rec = reactions[0]
    assert rec.reaction_smiles == "[CH3].[OH-]>[Na+]>[CH3OH]"
    assert rec.reactants == ["[CH3]", "[OH-]"]
    assert rec.reagents == ["[Na+]"]
    assert rec.products == ["[CH3OH]"]
    assert rec.atom_mapping == "[CH3:1].[OH-:2]>[Na+]>[CH3OH:1]"


def test_load_json_strip_atom_mapping(tmp_path):
    payload = [{"rxn": "[CH3:1].[OH-:2]>>[CH3OH:1]", "tag": "mapped"}]
    path = tmp_path / "rxn_map.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    def mapper(entry):
        return ReactionRecord(
            reaction_smiles=entry["rxn"],
            extra_metadata={"tag": entry["tag"]},
            source="json",
        )

    reactions = loader.load_json(str(path), mapper, strip_atom_mapping=True)

    assert len(reactions) == 1
    rec = reactions[0]
    assert rec.reaction_smiles == "[CH3].[OH-]>>[CH3OH]"
    assert rec.reactants == ["[CH3]", "[OH-]"]
    assert rec.reagents == []
    assert rec.products == ["[CH3OH]"]
    assert rec.atom_mapping == "[CH3:1].[OH-:2]>>[CH3OH:1]"
    assert rec.extra_metadata["tag"] == "mapped"
    assert rec.source == "json"


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
