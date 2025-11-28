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
    
    def show(
        self,
        size: tuple[int, int] = (2400, 600),
        jupyter: bool = True,
        with_meta: bool = True,
        show_atom_map_numbers: bool = False,
    ):
        from rdkit.Chem import Draw, rdChemReactions

        rxn_smiles = (self.raw)

        rxn = rdChemReactions.ReactionFromSmarts(
            rxn_smiles,
            useSmiles=True
        )

        if not show_atom_map_numbers:
            # Remove atom-map numbers to avoid numbered atom labels in depictions.
            for mol in (
                list(rxn.GetReactants())
                + list(rxn.GetProducts())
                + list(rxn.GetAgents())
            ):
                for atom in mol.GetAtoms():
                    if atom.HasProp("molAtomMapNumber"):
                        atom.ClearProp("molAtomMapNumber")

        img = Draw.ReactionToImage(
            rxn,
            subImgSize=(size[0] // 3, size[1]) 
        )

        if jupyter:
            try:
                from IPython.display import display, HTML
            except ImportError:  # pragma: no cover - optional dependency
                jupyter = False
            else:
                display(img)
                if with_meta and self.meta:
                    rows = "".join(
                        f"<tr><th>{k}</th><td>{v}</td></tr>"
                        for k, v in self.meta.items()
                    )
                    display(HTML(
                        "<table border='1' style='border-collapse:collapse;'>"
                        f"{rows}</table>"
                    ))

        if not jupyter:
            try:
                img.show()
            except Exception as exc:  # pragma: no cover - environment specific
                print(
                    "Unable to open the image viewer."
                    " Consider running in a Jupyter environment."
                )
                print(f"Original error: {exc}")
            if with_meta and self.meta:
                print("Metadata:")
                for k, v in self.meta.items():
                    print(f"  {k}: {v}")

        return img 
    

@dataclass
class ElementFilterRule:
    reactantElements: List[str]
    reagentElements: List[str]
    productElements: List[str]
