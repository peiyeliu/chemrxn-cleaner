# chemrxn_cleaner/loader.py

from __future__ import annotations
from typing import Callable, List, Dict, Any, Optional, Tuple

from ord_schema.message_helpers import (
    load_message,
    get_reaction_smiles,
)
from ord_schema.proto import dataset_pb2


def load_uspto_rsmi(
    path: str,
    keep_meta: bool = False,
) -> List[Tuple[str, Dict[str, Any]]]:
    """
    Load a USPTO .rsmi file.

    Typical format:
        reaction_smiles[\\t extra_field1 \\t extra_field2 ...]

    Args:
        path: Path to .rsmi file.
        keep_meta:
            - False: return List[(reaction_smiles, {})]
            - True:  return List[(reaction_smiles, meta_dict)]

    Returns:
        List[(str, dict)]
    """
    reactions_with_meta: List[Tuple[str, Dict[str, Any]]] = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            rxn_smiles = parts[0].strip()
            meta: Dict[str, Any] = {}
            if keep_meta and len(parts) > 1:
                meta["fields"] = parts[1:]
            reactions_with_meta.append((rxn_smiles, meta))

    return reactions_with_meta



def load_ord_pb_reaction_smiles(
    path: str,
    generate_if_missing: bool = True,
    allow_incomplete: bool = True,
    canonical: bool = True,
    meta_extractor: Optional[Callable[[dataset_pb2.Reaction], Dict[str, Any]]] = None,
) -> List[Tuple[str, Dict[str, Any]]]:
    """
    Load an ORD *.pb or *.pb.gz file and extract reaction SMILES.

    Args:
        path: Path to the ORD dataset file (e.g. 'data-00001.pb.gz').
        generate_if_missing:
            If True, generate reaction SMILES from inputs/outcomes when missing.
        allow_incomplete:
            If True, allow reactions with incomplete information when generating SMILES.
        canonical:
            If True, return canonical reaction SMILES.

    Returns:
        A list of (reaction_smiles, meta_dict) tuples.
    """
    dataset = load_message(path, dataset_pb2.Dataset)
    rxn_smiles_list: List[Tuple[str, Dict[str, Any]]] = []

    for idx, rxn in enumerate(dataset.reactions):
        smi: Optional[str] = get_reaction_smiles(
            message=rxn,
            generate_if_missing=generate_if_missing,
            allow_incomplete=allow_incomplete,
            canonical=canonical,
        )
        if smi:
            meta: Dict[str, Any] = {}
            if rxn.reaction_id:
                meta["reaction_id"] = rxn.reaction_id
            meta["reaction_index"] = idx

            if meta_extractor is not None:
                extra_meta = meta_extractor(rxn)
                if extra_meta:
                    meta.update(extra_meta)

            rxn_smiles_list.append((smi, meta))

    return rxn_smiles_list
