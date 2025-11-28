# chemrxn_cleaner/__init__.py


__version__ = "0.0.4"

import sys

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
from .io import (
    register_input_format,
    get_input_format,
    load_reactions,
    export_reaction_records_to_json,
    export_reaction_records_to_csv,
)
from .io import loader, loader_registry
from .extractor import ord_procedure_yields_meta

from .utils import similarity_filter
from .utils import similarity as _similarity

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
    "register_input_format",
    "get_input_format",
    "load_reactions",
    "export_reaction_records_to_json",
    "export_reaction_records_to_csv",
    # extractor
    "ord_procedure_yields_meta",
]

# Backward-compatible aliases for relocated submodules
sys.modules[__name__ + ".loader"] = loader
sys.modules[__name__ + ".loader_registry"] = loader_registry
sys.modules[__name__ + ".similarity"] = _similarity
