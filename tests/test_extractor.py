from chemrxn_cleaner import extractor


def test_ord_procedure_yields_meta(monkeypatch):
    class DummyProduct:
        def __init__(self, smiles: str, yield_value: float | None):
            self.smiles = smiles
            self.yield_value = yield_value

    class DummyOutcome:
        def __init__(self, products):
            self.products = products

    class DummyReaction:
        def __init__(self, outcomes, flat):
            self.outcomes = outcomes
            self.flat = flat

    reaction = DummyReaction(
        outcomes=[
            DummyOutcome(
                [
                    DummyProduct("P1", 90.0),
                    DummyProduct("P2", None),
                ]
            )
        ],
        flat={
            "setup.temperature": 298,
            "conditions.time": 5,
            "other.field": "skip",
        },
    )

    def fake_smiles_from_compound(product, canonical=True):
        return product.smiles

    def fake_get_product_yield(product, as_measurement=False):
        return product.yield_value

    def fake_message_to_row(rxn):
        return rxn.flat

    monkeypatch.setattr(extractor, "smiles_from_compound", fake_smiles_from_compound)
    monkeypatch.setattr(extractor, "get_product_yield", fake_get_product_yield)
    monkeypatch.setattr(extractor, "message_to_row", fake_message_to_row)

    meta = extractor.ord_procedure_yields_meta(reaction)

    assert meta == {
        "procedure": {
            "setup.temperature": 298,
            "conditions.time": 5,
        },
        "yields": [
            {"product_smiles": "P1", "yield_percent": 90.0},
        ],
    }
