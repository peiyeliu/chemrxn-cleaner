import pytest

from chemrxn_cleaner.cleaner import (
    clean_and_canonicalize,
    clean_reactions,
    clean_reactions_with_report,
)
from chemrxn_cleaner.parser import parse_reaction_smiles
from chemrxn_cleaner.types import ReactionRecord


def test_clean_reactions_attaches_metadata():
    rxns = [
        parse_reaction_smiles("C=CCBr>>C=CCI"),
        parse_reaction_smiles("CC(=O)O.OCC>[H+].[Cl-].OCC>CC(=O)OCC"),
        ReactionRecord(reaction_smiles="invalid-format"),
        parse_reaction_smiles("CC(=O)Cl.NH3>>CC(=O)NH2"),
    ]
    rxns[0].extra_metadata["source"] = "test"

    cleaned = clean_reactions(rxn_smiles_list=rxns, filters=[])

    assert len(cleaned) == 3
    assert cleaned[0].reaction_smiles == "C=CCBr>>C=CCI"
    assert cleaned[0].extra_metadata == {"source": "test"}
    assert cleaned[1].reaction_smiles == "CC(=O)O.OCC>[H+].[Cl-].OCC>CC(=O)OCC"
    assert cleaned[1].extra_metadata == {}
    assert cleaned[2].reaction_smiles == "CC(=O)Cl.NH3>>CC(=O)NH2"
    assert cleaned[2].extra_metadata == {}


def test_clean_reactions_raises_when_drop_disabled():
    rxns = [
        ReactionRecord(
            reaction_smiles="invalid-format", extra_metadata={"source": "bad"}
        ),
    ]

    with pytest.raises(ValueError):
        clean_reactions(rxn_smiles_list=rxns, drop_failed_parse=False, filters=[])


def test_clean_and_canonicalize_preserves_meta():
    rxns = [
        parse_reaction_smiles("OCC>O>CCO"),
    ]
    rxns[0].extra_metadata = {"id": 1}

    result = clean_and_canonicalize(rxn_smiles_list=rxns, isomeric=False, filters=[])

    assert len(result) == 1
    record = result[0]
    assert record.extra_metadata == {"id": 1}
    assert record.reactants == ["CCO"]
    assert record.products == ["CCO"]


def test_clean_reactions_with_report_tracks_stats():
    rxns = [
        ReactionRecord(reaction_smiles="CCO>>CCO"),
        ReactionRecord(reaction_smiles="CCO>>"),
        ReactionRecord(reaction_smiles="invalid-format"),
    ]

    cleaned, stats = clean_reactions_with_report(rxn_smiles_list=rxns)

    assert len(cleaned) == 1
    assert stats.n_input == 3
    assert stats.n_output == 1
    assert stats.n_failed_parse == 1

    has_product_stats = stats.per_filter["has_product"]
    assert has_product_stats.applied == 2
    assert has_product_stats.passed == 1
    assert has_product_stats.failed == 1

    molecules_valid_stats = stats.per_filter["all_molecules_valid"]
    assert molecules_valid_stats.applied == 1
    assert molecules_valid_stats.passed == 1
