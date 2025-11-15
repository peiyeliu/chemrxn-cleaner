# chemrxn-cleaner

A lightweight toolkit for cleaning and standardizing organic reaction datasets for machine learning.

## Prerequisites

- Python 3.9+
- [RDKit](https://www.rdkit.org/) (installable from `rdkit-pypi` on PyPI)
- Basic scientific Python stack (`pandas`, `tqdm`) — installed automatically via `pip install chemrxn-cleaner`
- Optional: `ord-schema` for loading Open Reaction Database datasets (`pip install chemrxn-cleaner[ord]`)

## Installation

```bash
pip install chemrxn-cleaner
# or from source
pip install -e .
```

## Quick Start

```python
from chemrxn_cleaner.loader import load_uspto_rsmi
from chemrxn_cleaner import basic_cleaning_pipeline, reporting

# Load reaction SMILES + metadata tuples
rxns = load_uspto_rsmi("/path/to/file.rsmi", keep_meta=True)

# Clean and canonicalize reactions
cleaned = basic_cleaning_pipeline(rxns)

# Summarize the cleaning process
report = reporting.summarize_cleaning(rxns, cleaned)
report.pretty_print()
```

### Loading ORD datasets

```python
from chemrxn_cleaner.loader import extract_ord_reaction_smiles_procedure_yield

rxns = extract_ord_reaction_smiles_procedure_yield("/path/to/ord_dataset.pb.gz")
```

### Custom Filtering with Metadata

```python
from chemrxn_cleaner.cleaning import clean_reactions
from chemrxn_cleaner.filters import meta_filter

rxns = load_uspto_rsmi("/path/to/file.rsmi", keep_meta=True)

# Keep only reactions where custom metadata flag is set
filters = [meta_filter(lambda meta: meta.get("is_valid", True))]
cleaned = clean_reactions(rxns, filters=filters)
```
