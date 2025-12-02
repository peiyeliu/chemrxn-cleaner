# chemrxn_cleaner/__init__.py


__version__ = "0.1.0"

import sys

from .types import ReactionRecord
from .parser import parse_reaction_smiles, canonicalize_reaction
from .filters import (
    ReactionFilter,
    has_product,
    all_molecules_valid,
    max_smiles_length,
    element_filter,
    meta_filter,
    default_filters,
)
from .cleaner import (
    clean_reactions,
    clean_and_canonicalize,
    basic_cleaning_pipeline,
)
from .ml import ForwardReactionDataset, records_to_dataframe, train_valid_test_split
from .io import (
    register_input_format,
    get_input_format,
    load_reactions,
    export_reaction_records_to_json,
    export_reaction_records_to_csv,
)
from .io import loader, loader_registry
from . import cleaner as _cleaner
from . import reporter as _reporter

from . import parser as _parser
from .utils import similarity_filter
from .utils import similarity as _similarity

__all__ = [
    # types
    "ReactionRecord",
    # parser
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
    # cleaner
    "clean_reactions",
    "clean_and_canonicalize",
    "basic_cleaning_pipeline",
    # io
    "register_input_format",
    "get_input_format",
    "load_reactions",
    "export_reaction_records_to_json",
    "export_reaction_records_to_csv",
    # ml
    "ForwardReactionDataset",
    "records_to_dataframe",
    "train_valid_test_split",
]
