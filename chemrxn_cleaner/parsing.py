# chemrxn_cleaner/parsing.py

from typing import Dict, List


def parse_reaction_smiles(rxn_smiles: str) -> Dict[str, List[str]]:
    """
    Parse a reaction SMILES of the form:
        reactants>reagents>products

    Returns:
        dict with keys: reactants, reagents, products
    """
    parts = rxn_smiles.split(">")
    if len(parts) != 3:
        raise ValueError(f"Invalid reaction SMILES format: {rxn_smiles!r}")

    reactants, reagents, products = parts
    return {
        "reactants": reactants.split(".") if reactants else [],
        "reagents": reagents.split(".") if reagents else [],
        "products": products.split(".") if products else [],
    }
