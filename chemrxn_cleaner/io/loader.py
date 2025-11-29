# chemrxn_cleaner/io/loader.py

from __future__ import annotations
import csv
import json

from typing import Callable, List, Dict, Any, Optional, Tuple, Sequence

from ord_schema.message_helpers import (
    load_message,
    get_reaction_smiles,
    get_product_yield,
    message_to_row,
)
from ord_schema.proto import dataset_pb2
from ..types import ReactionRecord, YieldType
from ..parser import parse_reaction_smiles


def load_uspto(
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


def load_ord(
    path: str,
    generate_if_missing: bool = True,
    allow_incomplete: bool = True,
    canonical: bool = True,
    meta_extractor: Optional[Callable[[dataset_pb2.Reaction], Dict[str, Any]]] = None,
) -> List[ReactionRecord]:
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
        A list of ReactionRecord objects with metadata + basic conditions populated.
    """
    dataset = load_message(path, dataset_pb2.Dataset)
    rxn_records: List[ReactionRecord] = []

    def _coerce_float(value: Any) -> Optional[float]:
        if value is None:
            return None
        if isinstance(value, (int, float)):
            return float(value)
        nested = getattr(value, "value", None)
        if nested is not None:
            try:
                return float(nested)
            except Exception:
                return None
        try:
            return float(value)
        except Exception:
            return None

    def _extract_conditions(
        reaction: dataset_pb2.Reaction,
    ) -> Dict[str, Any]:
        try:
            flat = message_to_row(reaction)
        except Exception:
            flat = {}

        def first(*keys: str) -> Any:
            for key in keys:
                if key in flat:
                    return flat[key]
            return None

        conditions: Dict[str, Any] = {}
        temperature = _coerce_float(
            first("conditions.temperature", "setup.temperature")
        )
        if temperature is not None:
            conditions["temperature_c"] = temperature

        time_hours = _coerce_float(first("conditions.time", "setup.time"))
        if time_hours is not None:
            conditions["time_hours"] = time_hours

        pressure_bar = _coerce_float(
            first("conditions.pressure", "setup.pressure", "conditions.vessel.pressure")
        )
        if pressure_bar is not None:
            conditions["pressure_bar"] = pressure_bar

        ph = _coerce_float(first("conditions.ph", "setup.ph"))
        if ph is not None:
            conditions["ph"] = ph

        atmosphere = first("conditions.atmosphere", "setup.atmosphere")
        if atmosphere is not None:
            conditions["atmosphere"] = str(atmosphere)

        scale_mmol = _coerce_float(first("conditions.scale", "setup.scale"))
        if scale_mmol is not None:
            conditions["scale_mmol"] = scale_mmol

        return conditions

    def _extract_primary_yield(
        reaction: dataset_pb2.Reaction,
    ) -> Tuple[Optional[float], YieldType]:
        best_yield: Optional[float] = None
        for outcome in getattr(reaction, "outcomes", []):
            for product in getattr(outcome, "products", []):
                try:
                    y_val = get_product_yield(product, as_measurement=False)
                except Exception:
                    y_val = None
                y_float = _coerce_float(y_val)
                if y_float is None:
                    continue
                if best_yield is None or y_float > best_yield:
                    best_yield = y_float
        if best_yield is None:
            return None, YieldType.NONE
        return best_yield, YieldType.OTHER

    for idx, rxn in enumerate(dataset.reactions):
        smi: Optional[str] = get_reaction_smiles(
            message=rxn,
            generate_if_missing=generate_if_missing,
            allow_incomplete=allow_incomplete,
            canonical=canonical,
        )
        if smi:
            record = parse_reaction_smiles(smi, strict=False)
            record.reaction_smiles = smi
            record.reaction_id = getattr(rxn, "reaction_id", "") or ""
            record.source = "ord"
            record.source_ref = path
            record.extra_metadata["reaction_index"] = idx

            yield_value, yield_type = _extract_primary_yield(rxn)
            record.yield_value = yield_value
            record.yield_type = yield_type

            for key, value in _extract_conditions(rxn).items():
                setattr(record, key, value)

            if meta_extractor is not None:
                extra_meta = meta_extractor(rxn)
                if extra_meta:
                    record.extra_metadata.update(extra_meta)

            rxn_records.append(record)

    return rxn_records


def load_csv(
    path: str,
    *,
    reactant_columns: Sequence[str] = (),
    product_columns: Sequence[str] = (),
    reagent_columns: Optional[Sequence[str]] = None,
    meta_columns: Optional[Sequence[str]] = None,
    reaction_smiles_column: Optional[str] = None,
    delimiter: str = ",",
) -> List[Tuple[str, Dict[str, Any]]]:
    """
    Load reaction SMILES assembled from CSV columns.

    Args:
        path: Path to the CSV file.
        reactant_columns: Column names that contain reactant SMILES.
        product_columns: Column names that contain product SMILES.
        reagent_columns: Optional column names that contain reagent SMILES.
        reaction_smiles_column:
            Optional column containing full reaction SMILES. If provided,
            reactant/product/reagent columns are ignored.
        meta_columns:
            Which columns to include as metadata. If None, use all columns
            not listed in reactant/product/reagent or reaction_smiles columns.
        delimiter: CSV delimiter (default ',').

    Returns:
        List[(str, dict)] where each tuple is (reaction_smiles, meta_dict).

    Raises:
        ValueError if required columns are missing.
    """

    def _normalize_columns(
        name: str, cols: Optional[Sequence[str]], allow_empty: bool = False
    ) -> List[str]:
        normalized = [str(col) for col in cols] if cols else []
        if not normalized and not allow_empty:
            raise ValueError(f"{name} must contain at least one column name.")
        return normalized

    def _join_smiles(row: Dict[str, Any], columns: List[str]) -> str:
        values: List[str] = []
        for col in columns:
            cell = row.get(col, "")
            if cell is None:
                continue
            text = str(cell).strip()
            if text:
                values.append(text)
        return ".".join(values)

    combined_mode = reaction_smiles_column is not None
    reactant_cols = _normalize_columns(
        "reactant_columns", reactant_columns, allow_empty=combined_mode
    )
    product_cols = _normalize_columns(
        "product_columns", product_columns, allow_empty=combined_mode
    )
    reagent_cols = _normalize_columns(
        "reagent_columns", reagent_columns, allow_empty=True
    )

    reactions: List[Tuple[str, Dict[str, Any]]] = []

    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        header = reader.fieldnames or []

        if combined_mode:
            if reaction_smiles_column not in header:
                raise ValueError(
                    f"CSV is missing required column: {reaction_smiles_column}"
                )
            used_columns = {reaction_smiles_column}
        else:
            required = set(reactant_cols + product_cols + reagent_cols)
            missing_required = sorted(required - set(header))
            if missing_required:
                missing_list = ", ".join(missing_required)
                raise ValueError(f"CSV is missing required columns: {missing_list}")
            used_columns = required

        if meta_columns is None:
            meta_cols = [col for col in header if col not in used_columns]
        else:
            meta_cols = _normalize_columns("meta_columns", meta_columns)
            missing_meta = sorted(set(meta_cols) - set(header))
            if missing_meta:
                missing_list = ", ".join(missing_meta)
                raise ValueError(
                    f"CSV is missing requested metadata columns: {missing_list}"
                )

        for row in reader:
            if combined_mode:
                rxn_smiles = str(row.get(reaction_smiles_column, "") or "").strip()
                if not rxn_smiles:
                    continue
            else:
                reactants_block = _join_smiles(row, reactant_cols)
                reagents_block = _join_smiles(row, reagent_cols)
                products_block = _join_smiles(row, product_cols)
                if not (reactants_block or reagents_block or products_block):
                    continue  # skip empty rows
                rxn_smiles = f"{reactants_block}>{reagents_block}>{products_block}"

            meta = {col: row.get(col, "") for col in meta_cols}
            reactions.append((rxn_smiles, meta))

    return reactions


def load_json(
    path: str,
    mapper: Callable[[Any], Optional[Tuple[str, Dict[str, Any]]]],
) -> List[Tuple[str, Dict[str, Any]]]:
    """
    Load reaction SMILES from a JSON file using a user-provided mapper.

    The JSON payload is expected to be a list; each entry is passed to
    ``mapper`` which should return ``(reaction_smiles, meta_dict)`` or ``None``
    to skip the entry.
    """
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)

    if not isinstance(payload, list):
        raise ValueError("JSON payload must be a list.")

    reactions: List[Tuple[str, Dict[str, Any]]] = []
    for idx, item in enumerate(payload):
        mapped = mapper(item)
        if mapped is None:
            continue
        if not (isinstance(mapped, (tuple, list)) and len(mapped) == 2):
            raise ValueError(
                f"Mapper must return (reaction_smiles, meta_dict); got {mapped!r} at index {idx}"
            )
        rxn_smiles, meta = mapped
        if not isinstance(rxn_smiles, str) or not rxn_smiles.strip():
            raise ValueError(
                f"Mapper must return a non-empty reaction_smiles at index {idx}"
            )
        if meta is None:
            meta = {}
        if not isinstance(meta, dict):
            raise ValueError(
                f"Mapper must return dict metadata at index {idx}, got {type(meta)!r}"
            )
        reactions.append((rxn_smiles.strip(), meta))

    return reactions
