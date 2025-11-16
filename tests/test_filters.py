import pytest

from chemrxn_cleaner.filters import meta_filter, element_filter
from chemrxn_cleaner.types import ElementFilterRule, ReactionRecord


def _record_with_meta(meta):
    return ReactionRecord(
        raw="CC(=O)O.OCC>[H+].[Cl-].OCC>CC(=O)OCC", 
        reactants=["CC(=O)O.OCC"], 
        reagents=["[H+].[Cl-].OCC"], 
        products=["CC(=O)OCC"], 
        meta=meta
    )


def test_meta_filter_uses_predicate_result():
    predicate = lambda meta: meta.get("score", 0) > 0.5
    filt = meta_filter(predicate)

    assert filt(_record_with_meta({"score": 0.7})) is True
    assert filt(_record_with_meta({"score": 0.1})) is False
    assert filt(_record_with_meta(None)) is False


def test_element_filter():
    rxn = _record_with_meta({})

    empty_rule = element_filter()

    forbid_Cl = element_filter(
        allowList=ElementFilterRule([],[],[]),
        forbidList=ElementFilterRule([],["Cl"],[])
    )

    allow_CH = element_filter(
        allowList=ElementFilterRule(["C", "H", "O", "Cl", "F", "N"],[],[]),
        forbidList=ElementFilterRule([],[],[])
    )

    assert empty_rule(rxn) is True
    assert forbid_Cl(rxn) is False
    assert allow_CH(rxn) is True


def test_element_filter_rejects_invalid_symbols():
    with pytest.raises(ValueError, match="Invalid element symbol"):
        element_filter(
            allowList=ElementFilterRule(["Xx"], [], []),
            forbidList=None,
        )


