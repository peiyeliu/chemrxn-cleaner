# chemrxn_cleaner/types.py
from dataclasses import dataclass
from typing import List, Dict, Any


@dataclass
class ReactionRecord:
    """Container for a single reaction instance."""
    raw: str                # raw SMILES strings
    reactants: List[str]
    reagents: List[str]
    products: List[str]
    meta: Dict[str, Any] | None = None 
