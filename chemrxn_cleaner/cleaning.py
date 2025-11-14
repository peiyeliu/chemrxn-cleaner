# chemrxn_cleaner/cleaning.py

from dataclasses import dataclass
from typing import Callable, List

from rdkit import Chem

from .parsing import parse_reaction_smiles


@dataclass
class ReactionRecord:
    raw: str
    reactants: List[str]
    reagents: List[str]
    products: List[str]


FilterFunc = Callable[[ReactionRecord], bool]


def is_valid_smiles(smiles: str) -> bool:
    """Return True if SMILES can be parsed by RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None


def default_valid_molecule_filter(smiles: str) -> bool:
    """Basic sanity check on a single molecule."""
    if not smiles:
        return False
    if not is_valid_smiles(smiles):
        return False
    # 简单长度限制，太长的分子先过滤掉（后面可以调参/改规则）
    if len(smiles) > 256:
        return False
    return True


def has_product_filter(record: ReactionRecord) -> bool:
    """Filter out reactions without products."""
    return len(record.products) > 0


def all_molecules_valid_filter(record: ReactionRecord) -> bool:
    """Filter out reactions with invalid SMILES."""
    all_smiles = record.reactants + record.reagents + record.products
    return all(default_valid_molecule_filter(s) for s in all_smiles)


def clean_reactions(
    rxn_smiles_list: List[str],
    filters: List[FilterFunc] | None = None,
) -> List[ReactionRecord]:
    """Apply parsing + filters to a list of reaction SMILES."""
    filters = filters or [has_product_filter, all_molecules_valid_filter]

    cleaned: List[ReactionRecord] = []
    for rxn in rxn_smiles_list:
        try:
            parsed = parse_reaction_smiles(rxn)
            record = ReactionRecord(raw=rxn, **parsed)
        except Exception:
            # 无法解析的直接丢弃
            continue

        if all(f(record) for f in filters):
            cleaned.append(record)

    return cleaned
