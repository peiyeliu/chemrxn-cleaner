# chemrxn_cleaner/filters.py

from __future__ import annotations

from typing import Callable, Iterable, List, Set, Dict, Any

from rdkit import Chem

from .types import ReactionRecord

ReactionFilter = Callable[[ReactionRecord], bool]


# ---------------------- tool functions ---------------------- #

def _iter_all_smiles(record: ReactionRecord) -> Iterable[str]:
    """Yield all SMILES strings from a ReactionRecord."""
    yield from record.reactants
    yield from record.reagents
    yield from record.products


def _is_valid_smiles(smiles: str) -> bool:
    """Return True if SMILES can be parsed by RDKit."""
    if not smiles:
        return False
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None


# ---------------------- basic filter ---------------------- #

def has_product(record: ReactionRecord) -> bool:
    """
    keep reactions have at least one product
    """
    return len(record.products) > 0


def all_molecules_valid(record: ReactionRecord) -> bool:
    """
    
    """
    for s in _iter_all_smiles(record):
        if not _is_valid_smiles(s):
            return False
    return True


# ---------------------- factory filters ---------------------- #

def max_smiles_length(max_len: int) -> ReactionFilter:
    """
    Ignore the molecule if its SMILES exceeds max_len
    """

    def _filter(record: ReactionRecord) -> bool:
        for s in _iter_all_smiles(record):
            if len(s) > max_len:
                return False
        return True

    return _filter


def allowed_elements(whitelist: Set[str]) -> ReactionFilter:
    """
    Return a filter function that keeps only reactions where every atom in
    every molecule is in the given `whitelist` of element symbols.
    """

    def _filter(record: ReactionRecord) -> bool:
        for s in _iter_all_smiles(record):
            if not s:
                return False
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                return False
            for atom in mol.GetAtoms():
                if atom.GetSymbol() not in whitelist:
                    return False
        return True

    return _filter


# ---------------------- metadata filter ---------------------- #



def meta_filter(predicate: Callable[[Dict[str, Any]], bool]) -> ReactionFilter:
    """
    Return a ReactionFilter that evaluates the provided predicate against the
    record metadata. Reactions with no metadata default to an empty dict.
    """

    def _filter(record: ReactionRecord) -> bool:
        meta = record.meta or {}
        try:
            return predicate(meta)
        except Exception:
            # Treat predicate errors as a failed filter match
            return False

    return _filter


# ---------------------- default filter ---------------------- #

def default_filters() -> List[ReactionFilter]:
    """
    return a default list of filters
    1. the reaction should contain at least one product
    2. all reactants and products should be valid
    """
    return [
        has_product,
        all_molecules_valid,
    ]
