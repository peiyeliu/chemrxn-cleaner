import pytest

from chemrxn_cleaner.cleaner import (
    clean_and_canonicalize,
    clean_reactions,
    clean_reactions_with_report,
)
from chemrxn_cleaner.parser import parse_reaction_smiles
from chemrxn_cleaner.reporter import CleaningStats, FilterStats
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


def test_cleaning_stats_combines_parallel_runs():
    stats_a = CleaningStats(n_input=3, n_output=2, n_failed_parse=1)
    stats_a.per_filter["has_product"] = FilterStats(
        name="has_product", applied=3, passed=2, failed=1
    )
    stats_a.per_filter["all_molecules_valid"] = FilterStats(
        name="all_molecules_valid", applied=1, passed=1, failed=0
    )

    stats_b = CleaningStats(n_input=2, n_output=1, n_failed_parse=0)
    stats_b.per_filter["has_product"] = FilterStats(
        name="has_product", applied=2, passed=1, failed=1
    )
    stats_b.per_filter["other_filter"] = FilterStats(
        name="other_filter", applied=5, passed=4, failed=1
    )

    combined = stats_a + stats_b

    assert combined.n_input == 5
    assert combined.n_output == 3
    assert combined.n_failed_parse == 1

    has_product = combined.per_filter["has_product"]
    assert has_product.applied == 5
    assert has_product.passed == 3
    assert has_product.failed == 2

    all_molecules_valid = combined.per_filter["all_molecules_valid"]
    assert all_molecules_valid.applied == 1
    assert all_molecules_valid.passed == 1
    assert all_molecules_valid.failed == 0

    other_filter = combined.per_filter["other_filter"]
    assert other_filter.applied == 5
    assert other_filter.passed == 4
    assert other_filter.failed == 1

    # Ensure operands are unchanged by pure addition
    assert stats_a.per_filter["has_product"].applied == 3
    assert stats_b.per_filter["has_product"].applied == 2

    stats_a += stats_b
    assert stats_a.n_input == 5
    assert stats_a.per_filter["has_product"].applied == 5
    assert stats_b.per_filter["has_product"].applied == 2


def test_filter_stats_str_format():
    fstats = FilterStats(name="foo", applied=3, passed=2, failed=1)
    assert str(fstats) == "foo: applied=3, passed=2, failed=1"


def test_cleaning_stats_summary_includes_filters():
    stats = CleaningStats(n_input=5, n_output=3, n_failed_parse=1)
    stats.per_filter["b_filter"] = FilterStats(
        name="b_filter", applied=2, passed=1, failed=1
    )
    stats.per_filter["a_filter"] = FilterStats(
        name="a_filter", applied=4, passed=3, failed=1
    )

    summary = stats.summary()

    assert "Input: 5, output: 3, failed_parse: 1" in summary
    lines = summary.splitlines()
    assert lines[1] == "Per-filter:"
    # Sorted order of keys should make the 'a_filter' line come first
    assert lines[2] == "  a_filter: applied=4, passed=3, failed=1"
    assert lines[3] == "  b_filter: applied=2, passed=1, failed=1"
    assert str(stats) == summary
