# chemrxn_cleaner/io/loader.py

from __future__ import annotations

import csv
import json
import logging
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from ord_schema.message_helpers import (
    get_product_yield,
    get_reaction_smiles,
    load_message,
    message_to_row,
)
from ord_schema.proto import dataset_pb2, reaction_pb2

from ..parser import parse_reaction_smiles
from ..types import ReactionRecord, YieldType

logger = logging.getLogger(__name__)


def load_uspto(
    path: str,
    keep_meta: bool = False,
) -> List[ReactionRecord]:
    """
    Load a USPTO .rsmi file.

    Typical format:
        reaction_smiles[\\t extra_field1 \\t extra_field2 ...]

    Args:
        path: Path to .rsmi file.
        keep_meta:
            - False: return List[ReactionRecord]
            - True:  return List[ReactionRecord] with extra_metadata["fields"] set

    Returns:
        List[ReactionRecord]
    """
    reactions: List[ReactionRecord] = []

    logger.info("Loading USPTO .rsmi file from %s", path)
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            rxn_smiles = parts[0].strip()
            if not rxn_smiles:
                continue

            record = parse_reaction_smiles(rxn_smiles, strict=False)
            record.reaction_smiles = rxn_smiles
            record.source = "uspto"
            record.source_file_path = path

            if keep_meta and len(parts) > 1:
                record.extra_metadata["fields"] = parts[1:]

            reactions.append(record)

    logger.info("Loaded %d reactions from USPTO file %s", len(reactions), path)
    return reactions


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
    logger.info(
        "Loading ORD dataset from %s (generate_if_missing=%s, allow_incomplete=%s, "
        "canonical=%s)",
        path,
        generate_if_missing,
        allow_incomplete,
        canonical,
    )
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

    def _extract_source_ref(reaction: dataset_pb2.Reaction) -> Optional[str]:
        provenance = getattr(reaction, "provenance", None)
        if provenance is None:
            return None

        for attr in ("doi", "patent", "publication_url"):
            value = getattr(provenance, attr, None)
            if value:
                return str(value)
        return None

    def _extract_procedure(
        reaction: reaction_pb2.Reaction,
    ) -> Dict[str, Any]:
        flat: Dict[str, Any] = message_to_row(reaction)
        procedure_prefixes = (
            "setup.",
            "conditions.",
            "workups.",
            "workup.",
            "notes.",
            "observations.",
        )
        return {
            k: v
            for k, v in flat.items()
            if any(k.startswith(pref) for pref in procedure_prefixes)
        }

    logger.debug(
        "Dataset contains %d reactions (path=%s)", len(dataset.reactions), path
    )
    for idx, rxn in enumerate(dataset.reactions):
        try:
            smi: Optional[str] = get_reaction_smiles(
                message=rxn,
                generate_if_missing=generate_if_missing,
                allow_incomplete=allow_incomplete,
                canonical=canonical,
            )
        except Exception:
            logger.exception(
                "Failed to extract reaction SMILES for reaction %s at index %d",
                getattr(rxn, "reaction_id", "") or "<unknown>",
                idx,
            )
            raise
        if smi:
            record = parse_reaction_smiles(smi, strict=False)
            record.reaction_smiles = smi
            record.reaction_id = getattr(rxn, "reaction_id", "") or ""
            record.source = "ord"
            record.source_file_path = path
            record.source_ref = _extract_source_ref(rxn)
            record.extra_metadata["reaction_index"] = idx

            yield_value, yield_type = _extract_primary_yield(rxn)
            record.yield_value = yield_value
            record.yield_type = yield_type

            for key, value in _extract_conditions(rxn).items():
                setattr(record, key, value)

            record.procedure = _extract_procedure(rxn)

            if meta_extractor is not None:
                extra_meta = meta_extractor(rxn)
                if extra_meta:
                    record.extra_metadata.update(extra_meta)

            rxn_records.append(record)
        else:
            logger.debug(
                "Skipping reaction at index %d due to missing reaction SMILES.", idx
            )

    logger.info("Loaded %d reactions from ORD dataset %s", len(rxn_records), path)
    return rxn_records


def load_csv(
    path: str,
    *,
    reactant_columns: Sequence[str] = (),
    product_columns: Sequence[str] = (),
    reagent_columns: Optional[Sequence[str]] = None,
    reaction_smiles_column: Optional[str] = None,
    delimiter: str = ",",
    skip_lines: int = 0,
    mapper: Optional[
        Callable[[ReactionRecord, Dict[str, Any]], Optional[ReactionRecord]]
    ] = None,
) -> List[ReactionRecord]:
    """
    Load ReactionRecord objects assembled from CSV columns.

    Args:
        path: Path to the CSV file.
        reactant_columns: Column names that contain reactant SMILES.
        product_columns: Column names that contain product SMILES.
        reagent_columns: Optional column names that contain reagent SMILES.
        reaction_smiles_column:
            Optional column containing full reaction SMILES. If provided,
            reactant/product/reagent columns are ignored.
        delimiter: CSV delimiter (default ',').
        skip_lines: Number of initial lines to skip before reading the header.
        mapper:
            Optional callable that receives the base ReactionRecord (parsed from the
            assembled reaction SMILES) and the raw row dictionary. It should return
            a ReactionRecord (possibly modified) or None to skip the row.

    Returns:
        List[ReactionRecord].

    Raises:
        ValueError if required columns are missing.
    """

    if skip_lines < 0:
        raise ValueError("skip_lines must be non-negative.")

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

    reactions: List[ReactionRecord] = []

    logger.info(
        "Loading CSV reactions from %s (combined_mode=%s, delimiter=%r, skip_lines=%d)",
        path,
        combined_mode,
        delimiter,
        skip_lines,
    )
    with open(path, "r", encoding="utf-8", newline="") as f:
        for _ in range(skip_lines):
            next(f, None)
        reader = csv.DictReader(f, delimiter=delimiter)
        header = reader.fieldnames or []

        if combined_mode:
            if reaction_smiles_column not in header:
                raise ValueError(
                    f"CSV is missing required column: {reaction_smiles_column}"
                )
        else:
            required = set(reactant_cols + product_cols + reagent_cols)
            missing_required = sorted(required - set(header))
            if missing_required:
                missing_list = ", ".join(missing_required)
                raise ValueError(f"CSV is missing required columns: {missing_list}")
        for idx, row in enumerate(reader):
            if combined_mode:
                rxn_smiles = str(row.get(reaction_smiles_column, "") or "").strip()
                if not rxn_smiles:
                    logger.debug("Skipping row %d: empty reaction_smiles column", idx)
                    continue
            else:
                reactants_block = _join_smiles(row, reactant_cols)
                reagents_block = _join_smiles(row, reagent_cols)
                products_block = _join_smiles(row, product_cols)
                if not (reactants_block or reagents_block or products_block):
                    logger.debug("Skipping row %d: no SMILES data found", idx)
                    continue  # skip empty rows
                rxn_smiles = f"{reactants_block}>{reagents_block}>{products_block}"

            record = parse_reaction_smiles(rxn_smiles, strict=False)
            record.reaction_smiles = rxn_smiles
            if not record.source:
                record.source = "csv"

            record.source_file_path = path

            if mapper is not None:
                mapped = mapper(record, row)
                if mapped is None:
                    continue
                if not isinstance(mapped, ReactionRecord):
                    raise ValueError(
                        f"Mapper must return a ReactionRecord;"
                        f"got {type(mapped)!r} at row {idx}"
                    )
                record = mapped

            if (
                not isinstance(record.reaction_smiles, str)
                or not record.reaction_smiles.strip()
            ):
                raise ValueError(
                    f"Mapper must set a non-empty reaction_smiles at row {idx}"
                )

            reactions.append(record)

    logger.info("Loaded %d reactions from CSV %s", len(reactions), path)
    return reactions


def load_json(
    path: str,
    mapper: Callable[[Any], Optional[ReactionRecord]],
) -> List[ReactionRecord]:
    """
    Load ReactionRecord instances from a JSON file using a user-provided mapper.

    The JSON payload is expected to be a list; each entry is passed to
    ``mapper`` which should return a populated ``ReactionRecord`` or ``None``
    to skip the entry.
    """
    logger.info("Loading JSON reactions from %s", path)
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)

    if not isinstance(payload, list):
        raise ValueError("JSON payload must be a list.")

    reactions: List[ReactionRecord] = []
    for idx, item in enumerate(payload):
        record = mapper(item)
        if record is None:
            continue
        if not isinstance(record, ReactionRecord):
            raise ValueError(
                f"Mapper must return a ReactionRecord;"
                f"got {type(record)!r} at index {idx}"
            )
        if (
            not isinstance(record.reaction_smiles, str)
            or not record.reaction_smiles.strip()
        ):
            raise ValueError(
                f"Mapper must set a non-empty reaction_smiles at index {idx}"
            )
        record.reaction_smiles = record.reaction_smiles.strip()
        reactions.append(record)

    logger.info("Loaded %d reactions from JSON %s", len(reactions), path)
    return reactions
