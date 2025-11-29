# similarity_filter.py

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Iterable, List, Literal, Union

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

from ..filters import ReactionFilter
from ..types import ReactionRecord


@dataclass(frozen=True)
class FingerprintConfig:
    radius: int = 2
    n_bits: int = 2048


# ====== Utility helpers: build RDKit molecules and fingerprints from SMILES ======


@lru_cache(maxsize=50_000)
def _mol_from_smiles(smiles: str) -> Chem.Mol | None:
    """Cached molecule parsing to avoid rebuilding the same structure."""
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return mol


@lru_cache(maxsize=50_000)
def _fingerprint_from_smiles(
    smiles: str, cfg: FingerprintConfig
) -> DataStructs.ExplicitBitVect | None:
    """Generate a Morgan fingerprint from SMILES and cache the result."""
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
    """Return the SMILES strings for molecules corresponding to the given role."""
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
    """
    Return a ReactionFilter for use in clean_reactions(..., filters=[...]).

    :param query_smiles: Target structure(s), one or more SMILES strings
    :param role: Which molecular role(s) in the reaction to match; "any" checks all three
    :param threshold: Tanimoto similarity cutoff; >= threshold counts as a hit
    :param radius: Morgan fingerprint radius (2 -> ECFP4, 3 -> ECFP6)
    :param n_bits: Fingerprint bit vector length
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
        """
        For a single reaction:
        - Collect all molecule SMILES for the configured role
        - Compare Tanimoto similarity with each query
        - Return True if any similarity >= threshold
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

    return _filter
