# chemrxn_cleaner/filters.py

from __future__ import annotations

from typing import Callable, Iterable, List, Set, Dict, Any

from rdkit import Chem

from .types import ReactionRecord, ElementFilterRule

ReactionFilter = Callable[[ReactionRecord], bool]

_PERIODIC_TABLE = Chem.GetPeriodicTable()


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


def _normalize_element_list(elements: Iterable[str] | None) -> Set[str] | None:
    """Return a set of elements, or None when no restriction is specified."""
    if not elements:
        return None
    normalized: Set[str] = set()
    for element in elements:
        if element is None:
            raise ValueError("Element entries must be strings.")
        symbol = element.strip()
        if not symbol:
            raise ValueError("Element symbol cannot be empty.")
        try:
            atomic_number = _PERIODIC_TABLE.GetAtomicNumber(symbol)
        except Exception as exc:
            raise ValueError(f"Invalid element symbol: {element!r}") from exc
        normalized.add(_PERIODIC_TABLE.GetElementSymbol(atomic_number))
    return normalized


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


def element_filter(
    allowList: ElementFilterRule | None = None,
    forbidList: ElementFilterRule | None = None,
) -> ReactionFilter:
    """
    Build a ReactionFilter that enforces per-role element rules.

    For each reactant, reagent, and product molecule, the generated filter
    ensures: (1) the SMILES string is parseable by RDKit, (2) every atom's
    symbol appears in the corresponding `allowList` set (when provided), and (3)
    no atom's symbol appears in the matching `forbidList` set (when provided).
    The reaction is rejected as soon as any molecule violates those
    constraints. Omitting either rule disables that portion of the filter.
    """

    allowed_reactants = _normalize_element_list(
        allowList.reactantElements if allowList else None
    )
    allowed_reagents = _normalize_element_list(
        allowList.reagentElements if allowList else None
    )
    allowed_products = _normalize_element_list(
        allowList.productElements if allowList else None
    )

    forbidden_reactants = _normalize_element_list(
        forbidList.reactantElements if forbidList else None
    )
    forbidden_reagents = _normalize_element_list(
        forbidList.reagentElements if forbidList else None
    )
    forbidden_products = _normalize_element_list(
        forbidList.productElements if forbidList else None
    )

    def _check_molecules(
        smiles_list: Iterable[str],
        allowed: Set[str] | None,
        forbidden: Set[str] | None,
    ) -> bool:
        for smile in smiles_list:
            mol = Chem.MolFromSmiles(smile)
            if mol is None:
                return False
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                if allowed is not None and symbol not in allowed:
                    return False
                if forbidden is not None and symbol in forbidden:
                    return False
        return True

    def _filter(record: ReactionRecord) -> bool:
        return (
            _check_molecules(record.reactants, allowed_reactants, forbidden_reactants)
            and _check_molecules(record.reagents, allowed_reagents, forbidden_reagents)
            and _check_molecules(record.products, allowed_products, forbidden_products)
        )

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
