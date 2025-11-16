"""
Example usage of chemrxn-cleaner on a tiny in-memory dataset.

Run with:
    python examples/clean_data_example.py
"""

from __future__ import annotations
import json
from chemrxn_cleaner.loader import load_uspto_rsmi, load_ord_pb_reaction_smiles
from chemrxn_cleaner.extractor import ord_procedure_yields_meta
from chemrxn_cleaner.cleaning import clean_reactions
from chemrxn_cleaner.filters import default_filters, max_smiles_length, element_filter
from chemrxn_cleaner.reporting import summarize_cleaning
from chemrxn_cleaner.types import ElementFilterRule

# read uspto data from your local environment
rxn_uspto = load_uspto_rsmi("/home/pyl/datasets/uspto/1976_Sep2016_USPTOgrants_smiles.rsmi")
cleaned_rxn_uspto = clean_reactions(rxn_smiles_list=rxn_uspto, filters=[])
clearning_report_uspto = summarize_cleaning(raw_reactions=rxn_uspto, cleaned_reactions=cleaned_rxn_uspto)


# read ORD (Open Reaction Database) data from your local environment
rxn_ord = load_ord_pb_reaction_smiles("/home/pyl/datasets/ord-data/data/00/ord_dataset-00005539a1e04c809a9a78647bea649c.pb.gz",
                                      meta_extractor=ord_procedure_yields_meta)
cleaned_rxn_ord = clean_reactions(rxn_smiles_list=rxn_ord, filters=[
    max_smiles_length(100),
    element_filter(
        forbidList=ElementFilterRule(["Cl", "Br"], [], [])
    )
])
report_ord = summarize_cleaning(raw_reactions=rxn_ord, cleaned_reactions=cleaned_rxn_ord)

# export results
with open('output.txt', 'w') as f:
    json.dump([rec.to_dict() for rec in cleaned_rxn_ord], f, indent=4) 
