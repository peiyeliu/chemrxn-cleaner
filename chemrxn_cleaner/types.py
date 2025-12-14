# chemrxn_cleaner/types.py
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional


class YieldType(str, Enum):
    """Yield measurement types supported by the dataset."""

    NONE = "none"  # Not a real yield
    ISOLATED = "isolated"  # isolated yield after workup and purification.
    ASSAY = "assay"  # Assay / crude yield (e.g., NMR, UPLC)
    LCAP = "lcap"  # LC area percent
    CONVERSION = "conversion"  # Conversion-based “yield” (GC, NMR, etc.)
    OTHER = "other"  # Non-standard yield type


@dataclass
class ReactionRecord:
    """Container for a single reaction instance.

    Attributes:
        reaction_id: Unique identifier for the reaction within a dataset.
        source: Origin of the reaction data (e.g., ``"uspto"``, ``"ord"``).
        source_ref: Reference such as DOI or patent number.
        source_file_path: File path from which the record was loaded.
        reaction_smiles: Raw reaction SMILES string.
        reactants: Reactant SMILES strings.
        reagents: Reagent/agent SMILES strings.
        products: Product SMILES strings.
        atom_mapping: Optional atom-mapped reaction SMILES.
        reaction_class: Optional reaction class label.
        procedure: Experimental procedure or flattened metadata.
        temperature_c: Reaction temperature in Celsius.
        time_hours: Reaction time in hours.
        pressure_bar: Reaction pressure in bar.
        ph: Observed or targeted pH.
        solvents: List of solvent identifiers.
        catalysts: List of catalysts.
        bases: List of bases.
        additives: List of other additives.
        atmosphere: Reaction atmosphere description.
        scale_mmol: Reaction scale in millimoles.
        yield_value: Reported yield numeric value.
        yield_type: Category describing how ``yield_value`` was measured.
        success: Optional boolean flag for reaction success.
        selectivity: Selectivity metric value, if provided.
        selectivity_type: Descriptor for the selectivity metric.
        is_balanced: Whether the reaction is atom-balanced.
        sanity_check_passed: Indicates whether sanity checks passed.
        warnings: Free-form warnings accumulated during processing.
        split: Dataset split label (e.g., train/valid/test).
        extra_numeric: Additional numeric metadata keyed by name.
        extra_categorical: Additional categorical metadata keyed by name.
        extra_metadata: Unstructured metadata not captured elsewhere.
    """

    # ---- Required identifiers / core representation ----
    reaction_id: str = ""
    source: str = ""
    source_ref: Optional[str] = None
    source_file_path: Optional[str] = None
    reaction_smiles: str = ""  # raw SMILES strings
    reactants: List[str] = field(default_factory=list)
    reagents: List[str] = field(default_factory=list)
    products: List[str] = field(default_factory=list)

    atom_mapping: Optional[str] = None
    reaction_class: Optional[str] = None

    # ---- Experimental conditions ----
    temperature_c: Optional[float] = None
    time_hours: Optional[float] = None
    pressure_bar: Optional[float] = None
    ph: Optional[float] = None
    procedure: Dict[str, Any] = field(default_factory=dict)

    solvents: List[str] = field(default_factory=list)
    catalysts: List[str] = field(default_factory=list)
    bases: List[str] = field(default_factory=list)
    additives: List[str] = field(default_factory=list)
    atmosphere: Optional[str] = None
    scale_mmol: Optional[float] = None

    # ---- Outcomes (targets / labels) ----
    yield_value: Optional[float] = None
    yield_type: YieldType = YieldType.NONE
    success: Optional[bool] = None

    # ---- selectivity ----
    selectivity: Optional[float] = None
    selectivity_type: Optional[str] = None

    # ---- Data-quality / bookkeeping ----
    is_balanced: Optional[bool] = None
    sanity_check_passed: bool = True
    warnings: List[str] = field(default_factory=list)

    # Dataset management only – e.g. “train/valid/test/time_split”
    split: Optional[str] = None

    # ---- Extension hooks ----
    extra_numeric: Dict[str, float] = field(default_factory=dict)
    extra_categorical: Dict[str, str] = field(default_factory=dict)
    extra_metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Return a plain-Python representation that is JSON serializable.

        Returns:
            Dictionary containing shallow copies of all fields suitable for JSON
            serialization.
        """
        return {
            "reaction_id": self.reaction_id,
            "source": self.source,
            "source_ref": self.source_ref,
            "source_file_path": self.source_file_path,
            "reaction_smiles": self.reaction_smiles,
            "reactants": list(self.reactants),
            "reagents": list(self.reagents),
            "products": list(self.products),
            "atom_mapping": self.atom_mapping,
            "reaction_class": self.reaction_class,
            "procedure": dict(self.procedure),
            "temperature_c": self.temperature_c,
            "time_hours": self.time_hours,
            "pressure_bar": self.pressure_bar,
            "ph": self.ph,
            "solvents": list(self.solvents),
            "catalysts": list(self.catalysts),
            "bases": list(self.bases),
            "additives": list(self.additives),
            "atmosphere": self.atmosphere,
            "scale_mmol": self.scale_mmol,
            "yield_value": self.yield_value,
            "yield_type": self.yield_type.value if self.yield_type else None,
            "success": self.success,
            "selectivity": self.selectivity,
            "selectivity_type": self.selectivity_type,
            "is_balanced": self.is_balanced,
            "sanity_check_passed": self.sanity_check_passed,
            "warnings": list(self.warnings),
            "split": self.split,
            "extra_numeric": dict(self.extra_numeric),
            "extra_categorical": dict(self.extra_categorical),
            "extra_metadata": dict(self.extra_metadata or {}),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ReactionRecord":
        """Construct a ``ReactionRecord`` from a dictionary payload.

        Args:
            data: Mapping produced by ``to_dict`` or a similarly structured
                dictionary.

        Returns:
            A populated ``ReactionRecord`` instance.
        """
        yield_type_value = data.get("yield_type") or YieldType.NONE
        if isinstance(yield_type_value, str):
            try:
                yield_type_enum = YieldType(yield_type_value)
            except ValueError:
                yield_type_enum = YieldType.NONE
        else:
            yield_type_enum = yield_type_value

        return cls(
            reaction_id=data.get("reaction_id", ""),
            source=data.get("source", ""),
            source_ref=data.get("source_ref"),
            source_file_path=data.get("source_file_path"),
            reaction_smiles=data.get("reaction_smiles") or data.get("raw", ""),
            reactants=list(data.get("reactants", []) or []),
            reagents=list(data.get("reagents", []) or []),
            products=list(data.get("products", []) or []),
            atom_mapping=data.get("atom_mapping"),
            reaction_class=data.get("reaction_class"),
            procedure=dict(data.get("procedure", {}) or {}),
            temperature_c=data.get("temperature_c"),
            time_hours=data.get("time_hours"),
            pressure_bar=data.get("pressure_bar"),
            ph=data.get("ph"),
            solvents=list(data.get("solvents", []) or []),
            catalysts=list(data.get("catalysts", []) or []),
            bases=list(data.get("bases", []) or []),
            additives=list(data.get("additives", []) or []),
            atmosphere=data.get("atmosphere"),
            scale_mmol=data.get("scale_mmol"),
            yield_value=data.get("yield_value"),
            yield_type=yield_type_enum,
            success=data.get("success"),
            selectivity=data.get("selectivity"),
            selectivity_type=data.get("selectivity_type"),
            is_balanced=data.get("is_balanced"),
            sanity_check_passed=data.get("sanity_check_passed", True),
            warnings=list(data.get("warnings", []) or []),
            split=data.get("split"),
            extra_numeric=dict(data.get("extra_numeric", {}) or {}),
            extra_categorical=dict(data.get("extra_categorical", {}) or {}),
            extra_metadata=cls._merge_extra_metadata(data),
        )

    @classmethod
    def _merge_extra_metadata(cls, data: Dict[str, Any]) -> Dict[str, Any]:
        """Combine explicit extra metadata with any unknown top-level keys."""
        base_meta = data.get("extra_metadata") or {}
        if isinstance(base_meta, dict):
            extra_meta = dict(base_meta)
        else:
            try:
                extra_meta = dict(base_meta)
            except Exception:
                extra_meta = {}

        meta_alias = data.get("meta")
        if isinstance(meta_alias, dict):
            extra_meta.update(meta_alias)

        known_keys = set(cls.__dataclass_fields__.keys()) | {"meta"}
        for key, value in data.items():
            if key not in known_keys and key not in extra_meta:
                extra_meta[key] = value

        return extra_meta

    def __post_init__(self) -> None:
        """Ensure optional mapping fields default to empty dictionaries."""
        if self.extra_metadata is None:
            self.extra_metadata = {}
        if self.procedure is None:
            self.procedure = {}

    def show(
        self,
        size: tuple[int, int] = (2400, 600),
        jupyter: bool = True,
        with_meta: bool = True,
        show_atom_map_numbers: bool = False,
    ):
        """Visualize the reaction using RDKit.

        Args:
            size: Tuple of pixel dimensions for the rendered image.
            jupyter: When True, display inline using IPython display hooks.
            with_meta: When True, show ``extra_metadata`` alongside the image.
            show_atom_map_numbers: When True, keep atom-map numbers in the
                depiction; when False, they are stripped before drawing.

        """
        from rdkit.Chem import Draw, rdChemReactions

        rxn_smiles = self.reaction_smiles

        rxn = rdChemReactions.ReactionFromSmarts(rxn_smiles, useSmiles=True)

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

        img = Draw.ReactionToImage(rxn, subImgSize=(size[0] // 3, size[1]))

        if jupyter:
            try:
                from IPython.display import HTML, display
            except ImportError:  # pragma: no cover - optional dependency
                jupyter = False
            else:
                display(img)
                if with_meta and self.extra_metadata:
                    rows = "".join(
                        f"<tr><th>{k}</th><td>{v}</td></tr>"
                        for k, v in self.extra_metadata.items()
                    )
                    display(
                        HTML(
                            "<table border='1' style='border-collapse:collapse;'>"
                            f"{rows}</table>"
                        )
                    )

        if not jupyter:
            try:
                img.show()
            except Exception as exc:  # pragma: no cover - environment specific
                print(
                    "Unable to open the image viewer."
                    " Consider running in a Jupyter environment."
                )
                print(f"Original error: {exc}")
            if with_meta and self.extra_metadata:
                print("Metadata:")
                for k, v in self.extra_metadata.items():
                    print(f"  {k}: {v}")

        return img
