# chemrxn-cleaner

Helpers for parsing, cleaning, filtering, reporting, and exporting organic reaction datasets before ML or analytics workflows.

## Installation

Requires Python 3.9+ with RDKit available in the environment (platform-specific wheels are not bundled). Install the package directly:

```bash
pip install chemrxn-cleaner
```

## Quick start

```python
from chemrxn_cleaner import (
    clean_and_canonicalize,
    clean_reactions_with_report,
    default_filters,
    export_reaction_records,
    load_reactions,
)

raw = load_reactions("data/sample.rsmi", fmt="uspto", keep_meta=True)
filters = default_filters()

cleaned = clean_and_canonicalize(raw, filters=filters)
cleaned_with_report, stats = clean_reactions_with_report(raw, filters=filters)
print(f"Input: {stats.n_input}, output: {stats.n_output}, failed: {stats.n_failed_parse}")

export_reaction_records(cleaned, "cleaned.json", fmt="json")
```

## Loading reaction data

Use the format registry instead of hand-rolled parsers:

```python
from chemrxn_cleaner import load_reactions

uspto_rxns = load_reactions("data/uspto_sample.rsmi", fmt="uspto", keep_meta=True)

ord_rxns = load_reactions(
    "data/ord_dataset.pb.gz",
    fmt="ord",
    generate_if_missing=True,
    allow_incomplete=True,
    canonical=True,
)
```

Register your own loader for custom formats:

```python
from chemrxn_cleaner import register_input_format
from chemrxn_cleaner.types import ReactionRecord

def load_my_format(path: str):
    rec = ReactionRecord(reaction_smiles="A>B>C", source="myfmt")
    return [rec]

register_input_format("myfmt", load_my_format)
rxns = load_reactions("my_file.txt", fmt="myfmt")
```

## Cleaning and reporting

Filters are plain callables. Compose the built-ins or author your own:

```python
from chemrxn_cleaner import (
    clean_and_canonicalize,
    clean_reactions_with_report,
    default_filters,
    max_smiles_length,
)

filters = default_filters() + [max_smiles_length(250)]
cleaned, stats = clean_reactions_with_report(raw, filters=filters)
canonicalized = clean_and_canonicalize(cleaned, filters=filters, isomeric=True)
print(stats.per_filter["max_smiles_length"].failed)
```

## Exports and ML utilities

```python
from chemrxn_cleaner import (
    ForwardReactionDataset,
    export_reaction_records,
    records_to_dataframe,
    train_valid_test_split,
)

df = records_to_dataframe(cleaned)
export_reaction_records(cleaned, "cleaned.csv", fmt="csv")

train, valid, test = train_valid_test_split(cleaned, seed=123)
dataset = ForwardReactionDataset(train, use_agents=True)
```

The API reference page documents every helper in detail.
