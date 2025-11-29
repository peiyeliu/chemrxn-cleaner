import pytest

from chemrxn_cleaner.cleaning import clean_reactions, clean_and_canonicalize


def test_clean_reactions_attaches_metadata():
    rxns = [
        ("C=CCBr>>C=CCI", {"source": "test"}),
        ("CC(=O)O.OCC>[H+].[Cl-].OCC>CC(=O)OCC", {}),
        ("invalid-format", {"source": "bad"}),
        ("CC(=O)Cl.NH3>>CC(=O)NH2", None),
    ]

    cleaned = clean_reactions(rxn_smiles_list=rxns, filters=[])

    assert len(cleaned) == 3
    assert cleaned[0].raw == "C=CCBr>>C=CCI"
    assert cleaned[0].meta == {"source": "test"}
    assert cleaned[1].raw == "CC(=O)O.OCC>[H+].[Cl-].OCC>CC(=O)OCC"
    assert cleaned[1].meta == {}
    assert cleaned[2].raw == "CC(=O)Cl.NH3>>CC(=O)NH2"
    assert cleaned[2].meta == {}


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
    assert record.meta == {"id": 1}
    assert record.reactants == ["CCO"]
    assert record.products == ["CCO"]
