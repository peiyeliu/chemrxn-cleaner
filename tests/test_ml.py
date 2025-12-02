import pandas as pd
import pytest

from chemrxn_cleaner import (
    ForwardReactionDataset,
    ReactionRecord,
    records_to_dataframe,
    train_valid_test_split,
)

torch = pytest.importorskip("torch")


def test_records_to_dataframe_flattens_lists_and_preserves_meta():
    records = [
        ReactionRecord(
            reaction_id="r1",
            reactants=["A", "B"],
            reagents=["NaH"],
            products=["C"],
            solvents=["THF"],
            warnings=["warn1", "warn2"],
            extra_metadata={"note": "batch1"},
        ),
        ReactionRecord(
            reaction_id="r2",
            reactants=["D"],
            reagents=[],
            products=["E"],
            solvents=[],
            extra_metadata={"note": "batch2"},
        ),
    ]

    df = records_to_dataframe(records)

    assert isinstance(df, pd.DataFrame)
    assert set(df["reaction_id"]) == {"r1", "r2"}
    row_r1 = df.loc[df["reaction_id"] == "r1"].iloc[0]
    assert row_r1["reactants"] == "A | B"
    assert row_r1["reagents"] == "NaH"
    assert row_r1["products"] == "C"
    assert row_r1["warnings"] == "warn1 | warn2"
    assert row_r1["solvents"] == "THF"
    assert row_r1["extra_metadata"] == {"note": "batch1"}


def test_train_valid_test_split_distributes_records_and_is_deterministic():
    records = [ReactionRecord(reaction_id=f"r{i}") for i in range(10)]

    train, valid, test = train_valid_test_split(records, seed=123)

    assert len(train) == 8
    assert len(valid) == 1
    assert len(test) == 1
    all_ids = {r.reaction_id for r in train + valid + test}
    assert all_ids == {f"r{i}" for i in range(10)}
    assert set(r.reaction_id for r in train).isdisjoint(
        set(r.reaction_id for r in valid + test)
    )
    assert set(r.reaction_id for r in valid).isdisjoint(
        set(r.reaction_id for r in test)
    )


def test_forward_reaction_dataset_builds_examples_with_and_without_agents():
    record = ReactionRecord(
        reaction_id="rxn-1",
        temperature_c=25.0,
        time_hours=2.5,
        solvents=["MeOH"],
        catalysts=["Pd/C"],
        bases=["K2CO3"],
        additives=["NaCl"],
        source="unit-test",
    )
    # ForwardReactionDataset expects these SMILES collections
    record.reactant_smiles = ["A", "B"]
    record.reagent_smiles = ["C"]
    record.product_smiles = ["D"]

    ds_with_agents = ForwardReactionDataset([record], use_agents=True)
    example = ds_with_agents[0]

    assert len(ds_with_agents) == 1
    assert example["input_smiles"] == "A.B.C"
    assert example["target_smiles"] == "D"
    assert example["reaction_id"] == "rxn-1"
    assert example["meta"]["temperature_c"] == 25.0
    assert example["meta"]["time_hours"] == 2.5
    assert example["meta"]["solvents"] == ["MeOH"]
    assert example["meta"]["catalysts"] == ["Pd/C"]
    assert example["meta"]["bases"] == ["K2CO3"]
    assert example["meta"]["additives"] == ["NaCl"]
    assert example["meta"]["source"] == "unit-test"

    ds_without_agents = ForwardReactionDataset([record], use_agents=False)
    example_no_agents = ds_without_agents[0]
    assert example_no_agents["input_smiles"] == "A.B"
