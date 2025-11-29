from chemrxn_cleaner.reporter import summarize_cleaning
from chemrxn_cleaner.types import ReactionRecord


def test_summarize_cleaning_counts_and_stats():
    raw = [
        ("rxn1", {"source": "a"}),
        ("rxn2", {"source": "b"}),
        ("rxn3", {"source": "c"}),
    ]
    cleaned = [
        ReactionRecord(
            reaction_smiles="rxn1",
            reactants=["A"],
            reagents=["B"],
            products=["C1", "C2"],
            extra_metadata={"source": "a"},
        ),
        ReactionRecord(
            reaction_smiles="rxn2",
            reactants=["A1", "A2"],
            reagents=[],
            products=["P"],
            extra_metadata={"source": "b"},
        ),
    ]

    report = summarize_cleaning(raw_reactions=raw, cleaned_reactions=cleaned)

    assert report.total_before == 3
    assert report.total_after == 2
    assert report.n_dropped() == 1
    assert report.avg_n_reactants == (1 + 2) / 2
    assert report.median_n_products == 1.5
