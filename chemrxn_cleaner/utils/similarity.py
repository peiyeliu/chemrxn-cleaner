# similarity_filter.py

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Iterable, List, Literal, Union

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from ..filters import ReactionFilter
from ..types import ReactionRecord


@dataclass(frozen=True)
class FingerprintConfig:
    """Configuration for generating Morgan fingerprints.

    Attributes:
        radius: Fingerprint radius (e.g., 2 for ECFP4).
        n_bits: Length of the fingerprint bit vector.
    """

    radius: int = 2
    n_bits: int = 2048


# ====== Utility helpers: build RDKit molecules and fingerprints from SMILES ======


@lru_cache(maxsize=50_000)
def _mol_from_smiles(smiles: str) -> Chem.Mol | None:
    """Parse SMILES into an RDKit molecule with caching.

    Args:
        smiles: Molecule SMILES string.

    Returns:
        RDKit molecule or ``None`` when parsing fails.
    """
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return mol


@lru_cache(maxsize=50_000)
def _fingerprint_from_smiles(
    smiles: str, cfg: FingerprintConfig
) -> DataStructs.ExplicitBitVect | None:
    """Generate a Morgan fingerprint from SMILES and cache the result.

    Args:
        smiles: Molecule SMILES string.
        cfg: Fingerprint configuration.

    Returns:
        Morgan fingerprint bit vector or ``None`` if parsing fails.
    """
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol,
        radius=cfg.radius,
        nBits=cfg.n_bits,
    )
    return fp


def _iter_molecule_smiles(rxn: ReactionRecord, role: str) -> Iterable[str]:
    """Yield molecule SMILES strings for the requested role.

    Args:
        rxn: Reaction record to inspect.
        role: One of ``"reactant"``, ``"reagent"``, ``"product"``, or ``"any"``.

    Returns:
        Iterator over SMILES strings in the specified role(s).

    Raises:
        ValueError: If ``role`` is not recognized.
    """
    if role == "reactant":
        yield from getattr(rxn, "reactants", [])
    elif role == "reagent":
        yield from getattr(rxn, "reagents", [])
    elif role == "product":
        yield from getattr(rxn, "products", [])
    elif role == "any":
        yield from getattr(rxn, "reactants", [])
        yield from getattr(rxn, "reagents", [])
        yield from getattr(rxn, "products", [])
    else:
        raise ValueError(f"Unknown role: {role}")


# ====== Public factory: build filter callable ======


def similarity_filter(
    query_smiles: Union[str, List[str]],
    role: Literal["reactant", "reagent", "product", "any"] = "any",
    threshold: float = 0.7,
    radius: int = 2,
    n_bits: int = 2048,
) -> ReactionFilter:
    """Return a similarity-based ``ReactionFilter``.

    Args:
        query_smiles: Target structure(s) expressed as SMILES.
        role: Molecular role(s) to check; ``"any"`` inspects reactants,
            reagents, and products.
        threshold: Tanimoto similarity cutoff; scores >= ``threshold`` count as
            a hit.
        radius: Morgan fingerprint radius (2 for ECFP4, 3 for ECFP6).
        n_bits: Length of the fingerprint bit vector.

    Returns:
        ReactionFilter predicate that keeps reactions containing at least one
        molecule similar to the query.

    Raises:
        ValueError: If any query SMILES cannot be parsed.
    """

    cfg = FingerprintConfig(radius=radius, n_bits=n_bits)

    # Normalize into a list
    if isinstance(query_smiles, str):
        query_list = [query_smiles]
    else:
        query_list = list(query_smiles)

    # Pre-compute all query fingerprints
    query_fps: List[DataStructs.ExplicitBitVect] = []
    for q in query_list:
        fp = _fingerprint_from_smiles(q, cfg)
        if fp is None:
            raise ValueError(f"Invalid query SMILES: {q}")
        query_fps.append(fp)

    def _filter(rxn: ReactionRecord) -> bool:
        """Return True when any molecule exceeds the similarity threshold.

        Args:
            rxn: Reaction to compare against the query fingerprints.

        Returns:
            True when at least one molecule meets or exceeds ``threshold``.
        """
        # Adjust here if you need every query matched instead of any
        for smi in _iter_molecule_smiles(rxn, role):
            fp = _fingerprint_from_smiles(smi, cfg)
            if fp is None:
                continue

            # Compare with all queries; any sufficiently similar molecule counts
            for qfp in query_fps:
                sim = DataStructs.TanimotoSimilarity(fp, qfp)
                if sim >= threshold:
                    return True

        return False

    _filter.__name__ = "similarity_filter"
    return _filter
