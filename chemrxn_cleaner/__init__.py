# chemrxn_cleaner/__init__.py


__version__ = "0.0.4"

from .types import ReactionRecord
from .parsing import parse_reaction_smiles, canonicalize_reaction
from .filters import (
    ReactionFilter,
    has_product,
    all_molecules_valid,
    max_smiles_length,
    element_filter,
    meta_filter,
    default_filters,
)
from .cleaning import (
    clean_reactions,
    clean_and_canonicalize,
    basic_cleaning_pipeline,
)
from .loader import (
    load_uspto_rsmi,
    load_ord_pb_reaction_smiles,
)
from .extractor import ord_procedure_yields_meta

from .similarity import similarity_filter

from .io import (
    export_reaction_records_to_json,
    load_reaction_records_from_json,
    export_reaction_records_to_csv,
    load_reaction_records_from_csv,
)

__all__ = [
    # types
    "ReactionRecord",
    # parsing
    "parse_reaction_smiles",
    "canonicalize_reaction",
    # filters
    "ReactionFilter",
    "has_product",
    "all_molecules_valid",
    "max_smiles_length",
    "element_filter",
    "meta_filter",
    "default_filters",
    "similarity_filter",
    # cleaning
    "clean_reactions",
    "clean_and_canonicalize",
    "basic_cleaning_pipeline",
    # io
    "load_uspto_rsmi",
    "load_ord_pb_reaction_smiles",
    "export_reaction_records_to_json",
    "load_reaction_records_from_json",
    "export_reaction_records_to_csv",
    "load_reaction_records_from_csv",
    # extractor
    "ord_procedure_yields_meta",
]
