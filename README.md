# chemrxn-cleaner

A lightweight toolkit for loading, cleaning, and standardizing organic reaction datasets before machine-learning or analytics workflows.

## Prerequisites

- Python 3.9+

These dependencies are pulled in automatically when installing `chemrxn-cleaner`, with the exception of platform-specific RDKit wheels.

## Installation

```bash
pip install chemrxn-cleaner
```

If you are working from a clone of this repository, install in editable mode to develop locally:

```bash
pip install -e .
```

## How to Use the Package

ChemRxn-Cleaner aims to make cleaning reproducible. A typical workflow has five steps:

1. Load reaction SMILES and metadata from USPTO `.rsmi` files or ORD protocol bundles.
2. Parse the SMILES strings into structured `ReactionRecord` objects.
3. Apply built-in or custom filters (length, element, metadata predicates, etc.).
4. Canonicalize the surviving reactions to obtain consistent representations.
5. Summarize the cleaning results or export the cleaned reactions downstream.

### 1. Loading Reaction Data

```python
from chemrxn_cleaner.loader import load_uspto_rsmi, load_ord_pb_reaction_smiles

# USPTO .rsmi loader (metadata fields stored in meta["fields"])
uspto_rxns = load_uspto_rsmi("data/uspto_sample.rsmi", keep_meta=True)

# ORD dataset loader with optional metadata extraction
from chemrxn_cleaner.extractor import ord_procedure_yields_meta
ord_rxns = load_ord_pb_reaction_smiles(
    "data/ord_dataset.pb.gz",
    meta_extractor=ord_procedure_yields_meta,
)
```

Both loaders return lists of `(reaction_smiles, metadata_dict)` tuples which feed directly into the cleaning utilities.

### 2. Running the Built-in Cleaning Pipeline

```python
from chemrxn_cleaner import basic_cleaning_pipeline

cleaned_ord = basic_cleaning_pipeline(ord_rxns)
```

`basic_cleaning_pipeline` parses, filters, and canonicalizes every reaction using the default filter stack (`has_product`, `all_molecules_valid`). The result is a list of immutable `ReactionRecord` instances with `reactants`, `reagents`, `products`, and `meta` attributes.

### 3. Custom Filters and Canonicalization Options

The cleaning helpers are composable; you can control the filter order and canonicalization behavior explicitly:

```python
from chemrxn_cleaner.cleaning import clean_and_canonicalize
from chemrxn_cleaner.filters import (
    default_filters,
    max_smiles_length,
    element_filter,
    meta_filter,
)
from chemrxn_cleaner.types import ElementFilterRule

filters = default_filters() + [
    max_smiles_length(250),
    element_filter(
        forbidList=ElementFilterRule(
            reactantElements=[],
            reagentElements=[],
            productElements=["Cl", "Br"],
        )
    ),
    meta_filter(lambda meta: meta.get("procedure", {}).get("setup.atmosphere") == "N2"),
]

cleaned_custom = clean_and_canonicalize(
    ord_rxns,
    filters=filters,
    isomeric=False,      # drop stereochemistry if desired
    drop_failed_parse=True,
)
```

Filters are simple callables accepting a `ReactionRecord` and returning `True`/`False`, so you can author domain-specific predicates without touching the core pipeline.

### 4. Working with Metadata

Metadata travels with each `ReactionRecord` through the pipeline. This makes it easy to slice cleaned reactions based on the original dataset context or to serialize extra descriptors:

```python
subset = [
    rec for rec in cleaned_custom
    if rec.meta.get("yields") and rec.meta["yields"][0]["yield_percent"] > 80
]
```

Custom metadata extractors (see `chemrxn_cleaner/extractor.py`) can capture procedure notes, yields, or any other ORD fields you care about.

### 5. Reporting and Exporting

```python
from chemrxn_cleaner import reporting

report = reporting.summarize_cleaning(ord_rxns, cleaned_custom)
report.pretty_print()

# Export canonical reaction SMILES + metadata for downstream use
import json
with open("cleaned_ord.jsonl", "w", encoding="utf-8") as f:
    for rec in cleaned_custom:
        payload = {
            "reactants": rec.reactants,
            "reagents": rec.reagents,
            "products": rec.products,
            "meta": rec.meta,
        }
        f.write(json.dumps(payload) + "\n")
```

## Quick Start

```python
from chemrxn_cleaner.loader import load_uspto_rsmi
from chemrxn_cleaner import basic_cleaning_pipeline, reporting

rxns = load_uspto_rsmi("/path/to/file.rsmi", keep_meta=True)
cleaned = basic_cleaning_pipeline(rxns)

report = reporting.summarize_cleaning(rxns, cleaned)
report.pretty_print()
```

## End-to-End Example Script

The repository ships with `examples/clean_data_example.py`, a runnable walkthrough that ties everything together:

1. Load USPTO `.rsmi` files and ORD `.pb.gz` datasets.
2. Apply default + custom filters (length, element whitelist/blacklist, metadata predicates).
3. Generate cleaning reports.
4. Persist cleaned reactions to disk.

Run it with:

```bash
python examples/clean_data_example.py
```

Use this script as a template—swap in your file paths, tweak the filter stack, and tailor the export code to match the format your downstream tools expect.
