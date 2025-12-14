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
    strip_atom_mapping: bool = False,
) -> List[ReactionRecord]:
    """Load a USPTO ``.rsmi`` file.

    The expected format is ``reaction_smiles[\\t extra_field1 \\t extra_field2 ...]``.

    Args:
        path: Path to the ``.rsmi`` file.
        keep_meta: When True, store trailing tab-separated fields under
            ``extra_metadata['fields']``.
        strip_atom_mapping: When True, remove atom-map numbers from input
            reaction SMILES and store the original mapped string in
            ``atom_mapping``.

    Returns:
        Parsed reaction records derived from the file.
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

            record = parse_reaction_smiles(
                rxn_smiles, strict=False, strip_atom_mapping=strip_atom_mapping
            )
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
    strip_atom_mapping: bool = False,
) -> List[ReactionRecord]:
    """Load an ORD ``*.pb`` or ``*.pb.gz`` file and extract reaction SMILES.

    Args:
        path: Path to the ORD dataset file (e.g., ``data-00001.pb.gz``).
        generate_if_missing: Generate reaction SMILES from inputs/outcomes when
            missing.
        allow_incomplete: Permit reactions with incomplete information when
            generating SMILES.
        canonical: Return canonical reaction SMILES when True.
        meta_extractor: Optional callable that produces extra metadata per
            reaction.
        strip_atom_mapping: When True, remove atom-map numbers from input
            reaction SMILES and store the original mapped string in
            ``atom_mapping``.

    Returns:
        Reaction records populated with parsed SMILES, metadata, and basic
        experimental conditions.
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
        """Convert numeric-like objects to floats when possible.

        Args:
            value: Input value that may encode a numeric measurement.

        Returns:
            Parsed float value or ``None`` when conversion is not possible.
        """
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
        """Extract a subset of reaction conditions from a protobuf reaction.

        Args:
            reaction: ORD reaction message.

        Returns:
            Dictionary containing standardized condition fields.
        """
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
        """Find the highest available yield and its type for a reaction.

        Args:
            reaction: ORD reaction message.

        Returns:
            Tuple of ``(yield_value, yield_type)`` where ``yield_value`` is
            ``None`` when no yield is available.
        """
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
        """Return a preferred provenance reference such as DOI or patent.

        Args:
            reaction: ORD reaction message.

        Returns:
            First available reference string, or ``None`` when absent.
        """
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
        """Collect flattened procedure-related fields for a reaction.

        Args:
            reaction: ORD reaction message.

        Returns:
            Dictionary containing only procedure-related keys.
        """
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
            record = parse_reaction_smiles(
                smi, strict=False, strip_atom_mapping=strip_atom_mapping
            )
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
    strip_atom_mapping: bool = False,
) -> List[ReactionRecord]:
    """Load ReactionRecord objects assembled from CSV columns.

    Args:
        path: Path to the CSV file.
        reactant_columns: Column names that contain reactant SMILES.
        product_columns: Column names that contain product SMILES.
        reagent_columns: Optional column names that contain reagent SMILES.
        reaction_smiles_column: Column containing full reaction SMILES; when
            provided, component columns are ignored.
        delimiter: CSV delimiter character.
        skip_lines: Number of initial lines to skip before reading the header.
        mapper: Optional callable that receives the parsed ReactionRecord and
            raw row dictionary. It may return a modified record or ``None`` to
            skip the row.
        strip_atom_mapping: When True, remove atom-map numbers from input
            reaction SMILES and store the original mapped string in
            ``atom_mapping``.

    Returns:
        Parsed reaction records derived from the CSV rows.

    Raises:
        ValueError: If required columns are missing or the mapper returns
            invalid data.
    """

    if skip_lines < 0:
        raise ValueError("skip_lines must be non-negative.")

    def _normalize_columns(
        name: str, cols: Optional[Sequence[str]], allow_empty: bool = False
    ) -> List[str]:
        """Validate and coerce a sequence of column names to strings.

        Args:
            name: Human-readable label for the columns being normalized.
            cols: Column names to normalize.
            allow_empty: Whether an empty list is permitted.

        Returns:
            Normalized list of column name strings.

        Raises:
            ValueError: If ``allow_empty`` is False and no columns are provided.
        """
        normalized = [str(col) for col in cols] if cols else []
        if not normalized and not allow_empty:
            raise ValueError(f"{name} must contain at least one column name.")
        return normalized

    def _join_smiles(row: Dict[str, Any], columns: List[str]) -> str:
        """Join SMILES fragments from selected columns with dots.

        Args:
            row: CSV row dictionary.
            columns: Ordered column names to collect SMILES fragments from.

        Returns:
            Dot-joined SMILES string (may be empty).
        """
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

            record = parse_reaction_smiles(
                rxn_smiles, strict=False, strip_atom_mapping=strip_atom_mapping
            )
            if not strip_atom_mapping:
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
    *,
    strip_atom_mapping: bool = False,
) -> List[ReactionRecord]:
    """Load ReactionRecord instances from a JSON list payload.

    Args:
        path: Path to the JSON file containing a list payload.
        mapper: Callable that converts each list item into a ``ReactionRecord``
            or ``None`` to skip it.
        strip_atom_mapping: When True, remove atom-map numbers from mapper-
            returned reaction SMILES and store the original mapped string in
            ``atom_mapping``.

    Returns:
        Reaction records produced by the mapper.

    Raises:
        ValueError: If the payload is not a list, the mapper returns invalid
            data, or a record is missing ``reaction_smiles``.
    """
    logger.info("Loading JSON reactions from %s", path)
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)

    if not isinstance(payload, list):
        raise ValueError("JSON payload must be a list.")

    reactions: List[ReactionRecord] = []

    def _strip_record(record: ReactionRecord) -> ReactionRecord:
        if not strip_atom_mapping:
            return record
        parsed = parse_reaction_smiles(
            record.reaction_smiles, strict=False, strip_atom_mapping=True
        )
        record.reaction_smiles = parsed.reaction_smiles
        record.reactants = parsed.reactants
        record.reagents = parsed.reagents
        record.products = parsed.products
        record.atom_mapping = parsed.atom_mapping or record.atom_mapping
        return record

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
        record = _strip_record(record)
        reactions.append(record)

    logger.info("Loaded %d reactions from JSON %s", len(reactions), path)
    return reactions
