"""
Helper utilities that extract structured metadata from ORD reactions.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from ord_schema.message_helpers import (
    get_product_yield,
    smiles_from_compound,
    message_to_row,
)
from ord_schema.proto import reaction_pb2


def ord_procedure_yields_meta(
    reaction: reaction_pb2.Reaction,
) -> Dict[str, Any]:
    """
    Build a metadata dict containing procedure fields and product yields
    for a single ORD reaction message.
    """
    yields_info: List[Dict[str, Any]] = []
    for outcome in reaction.outcomes:
        for product in outcome.products:
            prod_smi: Optional[str] = None
            try:
                prod_smi = smiles_from_compound(product, canonical=True)
            except Exception:
                prod_smi = None

            y_val = get_product_yield(product, as_measurement=False)
            if y_val is None:
                continue

            yields_info.append(
                {
                    "product_smiles": prod_smi,
                    "yield_percent": float(y_val),
                }
            )

    flat: Dict[str, Any] = message_to_row(reaction)
    procedure_prefixes = (
        "setup.",
        "conditions.",
        "workups.",
        "workup.",
        "notes.",
        "observations.",
    )
    procedure: Dict[str, Any] = {
        k: v
        for k, v in flat.items()
        if any(k.startswith(pref) for pref in procedure_prefixes)
    }

    return {
        "procedure": procedure,
        "yields": yields_info,
    }
