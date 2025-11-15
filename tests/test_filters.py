from chemrxn_cleaner.filters import meta_filter
from chemrxn_cleaner.types import ReactionRecord


def _record_with_meta(meta):
    return ReactionRecord(raw="", reactants=[], reagents=[], products=[], meta=meta)


def test_meta_filter_uses_predicate_result():
    predicate = lambda meta: meta.get("score", 0) > 0.5
    filt = meta_filter(predicate)

    assert filt(_record_with_meta({"score": 0.7})) is True
    assert filt(_record_with_meta({"score": 0.1})) is False
    assert filt(_record_with_meta(None)) is False
