# chemrxn-cleaner

Lightweight helpers for parsing, cleaning, filtering, and exporting organic reaction datasets before ML or analytics workflows.

## Installation

- Python 3.9+
- RDKit is required (platform-specific wheels are not bundled).

```bash
pip install chemrxn-cleaner
```

Developing locally? Install in editable mode:

```bash
pip install -e .
```

## Quick start

```python
from chemrxn_cleaner import basic_cleaning_pipeline, reporter, parse_reaction_smiles
from chemrxn_cleaner.io.loader import load_uspto

raw = load_uspto("data/sample.rsmi", keep_meta=True)

cleaned = basic_cleaning_pipeline(raw)

summary = reporter.summarize_cleaning(raw_reactions=raw, cleaned_reactions=cleaned)
summary.pretty_print()
```

## Loading reaction data

```python
from chemrxn_cleaner.io.loader import load_uspto, load_csv, load_ord, load_json
from chemrxn_cleaner.parser import parse_reaction_smiles

# USPTO .rsmi loader (optional metadata fields stored in extra_metadata["fields"])
uspto_rxns = load_uspto("data/uspto_sample.rsmi", keep_meta=True)

# CSV loader: assemble reaction SMILES from column mappings
csv_rxns = load_csv(
    "data/reactions.csv",
    reactant_columns=["reactant_a", "reactant_b"],
    reagent_columns=["catalyst"],
    product_columns=["product"],
    mapper=lambda record, row: (
        record.extra_metadata.update({"temperature": row.get("temp_c")}) or record
    ),
)

# CSV loader with a pre-built reaction_smiles column
csv_rxns_prebuilt = load_csv(
    "data/reactions.csv",
    reaction_smiles_column="rxn_smiles",
    mapper=lambda record, row: record,
)

The `mapper` callable receives `(record, row)` and can set optional attributes or skip rows by returning `None`.

# JSON loader with a custom mapper per entry
def map_json_entry(item):
    rec = parse_reaction_smiles(f"{item['reactants']}>>{item['products']}", strict=False)
    rec.source = "json"
    rec.extra_metadata.update(item.get("meta", {}))
    return rec

json_rxns = load_json("data/reactions.json", mapper=map_json_entry)

# ORD dataset loader (returns populated ReactionRecord objects)
ord_rxns = load_ord(
    "data/ord_dataset.pb.gz",
)
```

`load_uspto`, `load_csv`, `load_json`, and `load_ord` return `ReactionRecord` objects. `load_ord` additionally populates `reaction_id`, yields, basic conditions (temperature, time, pressure, pH, atmosphere, scale), and `extra_metadata["reaction_index"]`.

You can also register custom loaders and call them through the registry:

```python
from chemrxn_cleaner.io import register_input_format, load_reactions

def load_my_format(path: str):
    return [("A>B>C", {"source": path})]

register_input_format("myfmt", load_my_format)
rxns = load_reactions("my_file.txt", fmt="myfmt")
```

## Cleaning and filters

```python
from chemrxn_cleaner.cleaner import clean_reactions, clean_and_canonicalize
from chemrxn_cleaner.filters import (
    default_filters,
    max_smiles_length,
    element_filter,
    meta_filter,
    ElementFilterRule,
)
from chemrxn_cleaner.utils import similarity_filter

filters = default_filters() + [
    max_smiles_length(250),
    element_filter(
        forbidList=ElementFilterRule([], ["Cl"], []),
    ),
    meta_filter(lambda meta: meta.get("source") == "trusted"),
    similarity_filter("c1ccccc1", role="reactant", threshold=0.6),
]

cleaned = clean_and_canonicalize(
    rxn_smiles_list=uspto_rxns,
    filters=filters,
    isomeric=True,
)
```

- `clean_reactions` accepts `ReactionRecord` objects (parses `reaction_smiles` if reactants/reagents/products are empty) and applies filters.
- `clean_and_canonicalize` also canonicalizes every molecule; `basic_cleaning_pipeline` is the default stack (`has_product`, `all_molecules_valid`, strict parsing, isomeric SMILES).
- Filters are callables returning `True`/`False`; author your own to encode domain rules.

## Reporting and exporting

```python
from chemrxn_cleaner import reporter
from chemrxn_cleaner.io import export_reaction_records_to_json, export_reaction_records_to_csv

report = reporter.summarize_cleaning(raw, cleaned)
report.pretty_print()

export_reaction_records_to_json(cleaned, "cleaned.json")
export_reaction_records_to_csv(cleaned, "cleaned.csv")
```

## Working with `ReactionRecord`

`ReactionRecord` holds the parsed reaction (`reactants`, `reagents`, `products`, `reaction_smiles`) plus optional metadata like yields, conditions, and arbitrary `extra_metadata`. Use `to_dict()`/`from_dict()` for serialization, and `show()` to render a reaction image when RDKit visualization is available.

## Examples

An interactive walkthrough lives at `examples/example.ipynb`. Open it in Jupyter, swap in your own file paths, and adapt the filter stack for your dataset.
