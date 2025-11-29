import pytest

from chemrxn_cleaner.cleaner import clean_reactions, clean_and_canonicalize


def test_clean_reactions_attaches_metadata():
    rxns = [
        ("C=CCBr>>C=CCI", {"source": "test"}),
        ("CC(=O)O.OCC>[H+].[Cl-].OCC>CC(=O)OCC", {}),
        ("invalid-format", {"source": "bad"}),
        ("CC(=O)Cl.NH3>>CC(=O)NH2", None),
    ]

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
        ("invalid-format", {"source": "bad"}),
    ]

    with pytest.raises(ValueError):
        clean_reactions(rxn_smiles_list=rxns, drop_failed_parse=False, filters=[])


def test_clean_and_canonicalize_preserves_meta():
    rxns = [
        ("OCC>O>CCO", {"id": 1}),
    ]

    result = clean_and_canonicalize(rxn_smiles_list=rxns, isomeric=False, filters=[])

    assert len(result) == 1
    record = result[0]
    assert record.extra_metadata == {"id": 1}
    assert record.reactants == ["CCO"]
    assert record.products == ["CCO"]
