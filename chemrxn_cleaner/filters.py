# chemrxn_cleaner/filters.py

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, List, Set

from rdkit import Chem

from .types import ReactionRecord


@dataclass
class ElementFilterRule:
    """Allowed or forbidden element symbols by reaction role.

    Attributes:
        reactantElements: Allowed/forbidden element symbols for reactants.
        reagentElements: Allowed/forbidden element symbols for reagents.
        productElements: Allowed/forbidden element symbols for products.
    """

    reactantElements: List[str]
    reagentElements: List[str]
    productElements: List[str]


ReactionFilter = Callable[[ReactionRecord], bool]

_PERIODIC_TABLE = Chem.GetPeriodicTable()


# ---------------------- tool functions ---------------------- #


def _iter_all_smiles(record: ReactionRecord) -> Iterable[str]:
    """Yield all SMILES strings from a reaction record in canonical order.

    Args:
        record: Reaction to iterate through.

    Returns:
        Iterator over reactant, reagent, then product SMILES.
    """
    yield from record.reactants
    yield from record.reagents
    yield from record.products


def _is_valid_smiles(smiles: str) -> bool:
    """Check whether RDKit can parse the provided SMILES.

    Args:
        smiles: Candidate SMILES string.

    Returns:
        True when a molecule can be built from ``smiles``; otherwise False.
    """
    if not smiles:
        return False
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None


def _normalize_element_list(elements: Iterable[str] | None) -> Set[str] | None:
    """Normalize element symbols and return a unique set.

    Args:
        elements: Iterable of raw element symbols; ``None`` means no constraint.

    Returns:
        A set of cleaned element symbols, or ``None`` when ``elements`` is falsy.

    Raises:
        ValueError: If an element entry is empty or not recognized.
    """
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
    """Check whether the reaction contains at least one product molecule.

    Args:
        record: Reaction to evaluate.

    Returns:
        True when at least one product SMILES is present.
    """
    return len(record.products) > 0


def all_molecules_valid(record: ReactionRecord) -> bool:
    """Verify that all reaction molecules are valid SMILES.

    Args:
        record: Reaction to evaluate.

    Returns:
        True when every reactant, reagent, and product SMILES parses in RDKit.
    """
    for s in _iter_all_smiles(record):
        if not _is_valid_smiles(s):
            return False
    return True


# ---------------------- factory filters ---------------------- #


def max_smiles_length(max_len: int) -> ReactionFilter:
    """Build a filter that enforces a maximum SMILES length.

    Args:
        max_len: Maximum allowed character length for any SMILES string.

    Returns:
        ReactionFilter predicate that rejects reactions containing longer SMILES.
    """

    def _filter(record: ReactionRecord) -> bool:
        """Return True when all SMILES strings are within the length limit.

        Args:
            record: Reaction to evaluate.

        Returns:
            True if every SMILES length is less than or equal to ``max_len``.
        """
        for s in _iter_all_smiles(record):
            if len(s) > max_len:
                return False
        return True

    _filter.__name__ = f"max_smiles_length({max_len})"
    return _filter


def element_filter(
    allowList: ElementFilterRule | None = None,
    forbidList: ElementFilterRule | None = None,
) -> ReactionFilter:
    """Build a ReactionFilter enforcing per-role element rules.

    For each reactant, reagent, and product molecule, the generated filter:
    1) parses the SMILES with RDKit, 2) checks membership in ``allowList`` when
    provided, and 3) rejects any atom that appears in ``forbidList``.
    Constraints are evaluated per role; missing rules skip that check.

    Args:
        allowList: Allowed element symbols for each role. ``None`` disables
            allow-list filtering.
        forbidList: Forbidden element symbols for each role. ``None`` disables
            forbid-list filtering.

    Returns:
        ReactionFilter predicate implementing the configured constraints.

    Raises:
        ValueError: If element symbols are empty or invalid.
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
        """Validate molecules against allowed/forbidden element sets.

        Args:
            smiles_list: SMILES strings to validate.
            allowed: Elements permitted for the given role, or ``None`` to skip.
            forbidden: Elements to reject for the given role, or ``None`` to skip.

        Returns:
            True if all molecules satisfy both element constraints.
        """
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

    _filter.__name__ = "element_filter"
    return _filter


# ---------------------- metadata filter ---------------------- #


def meta_filter(predicate: Callable[[Dict[str, Any]], bool]) -> ReactionFilter:
    """Create a filter that tests metadata via a user predicate.

    Args:
        predicate: Callable that receives ``extra_metadata`` and returns True to
            keep the reaction or False to drop it.

    Returns:
        ReactionFilter wrapping the predicate.
    """

    def _filter(record: ReactionRecord) -> bool:
        """Return True when the predicate accepts the record metadata.

        Args:
            record: Reaction to evaluate.

        Returns:
            Result of the predicate; False when the predicate raises.
        """
        meta = record.extra_metadata or {}
        try:
            return predicate(meta)
        except Exception:
            # Treat predicate errors as a failed filter match
            return False

    pred_name = getattr(predicate, "__name__", predicate.__class__.__name__)
    _filter.__name__ = f"meta_filter({pred_name})"
    return _filter


# ---------------------- default filter ---------------------- #


def default_filters() -> List[ReactionFilter]:
    """Return the default list of reaction filters.

    The defaults ensure at least one product is present and that all SMILES
    strings are valid RDKit molecules.

    Returns:
        List of basic ReactionFilter callables.
    """
    return [
        has_product,
        all_molecules_valid,
    ]
