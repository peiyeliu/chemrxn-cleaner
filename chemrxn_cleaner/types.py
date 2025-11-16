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

    def to_dict(self) -> Dict[str, Any]:
        """
        Return a plain-Python representation that is JSON serializable.
        """
        return {
            "raw": self.raw,
            "reactants": list(self.reactants),
            "reagents": list(self.reagents),
            "products": list(self.products),
            "meta": self.meta,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ReactionRecord":
        """
        Construct a ReactionRecord from a dictionary structure (inverse of to_dict()).
        """
        return cls(
            raw=data.get("raw", ""),
            reactants=list(data.get("reactants", []) or []),
            reagents=list(data.get("reagents", []) or []),
            products=list(data.get("products", []) or []),
            meta=data.get("meta"),
        )
    

@dataclass
class ElementFilterRule:
    reactantElements: List[str]
    reagentElements: List[str]
    productElements: List[str]
