---
title: 'chemrxn-cleaner: A Lightweight Toolkit for Parsing, Cleaning, and Preparing Organic Reaction Data for Machine Learning'
tags:
  - chemistry
  - cheminformatics
  - reaction data
  - data cleaning
  - machine learning
  - Python
authors:
  - name: Your Name
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Independent Researcher, California, USA
    index: 1
date: 2025-XX-XX
bibliography: paper.bib
---

# Summary

Machine learning models for reaction prediction, yield estimation, and synthetic planning typically require large, diverse, and standardized reaction datasets. However, widely used data sources such as the USPTO reaction corpus, the Open Reaction Database (ORD), and custom laboratory exports vary significantly in structure, metadata availability, and quality. Before these datasets can be used for modeling or analytics, practitioners must address challenges such as inconsistent SMILES representations, missing role assignments, noisy or malformed entries, non-standard metadata formats, and heterogeneous file structures.

**chemrxn-cleaner** is a lightweight, extensible Python package that provides uniform tools for loading, parsing, cleaning, filtering, reporting, and exporting reaction data. It enables researchers to rapidly transform raw reaction datasets into standardized ML-ready representations with minimal boilerplate, while also supporting custom formats and user-defined domain rules. The package aims to bridge a critical but often under-documented part of the reaction-ML workflow: *data preparation*.

# Statement of Need

Reaction-focused ML research depends heavily on data quality. Even minor inconsistencies in reactant–reagent separation, SMILES formatting, or metadata propagation can degrade downstream model performance. Existing cheminformatics toolkits (e.g., RDKit) provide reaction parsing primitives but do not offer a complete, opinionated pipeline for dataset-level ingestion and cleaning. Meanwhile, large datasets such as ORD or USPTO require specialized logic to interpret their schema or metadata.

chemrxn-cleaner addresses these gaps by offering:

- A **unified interface** to load heterogeneous data sources (USPTO `.rsmi`, ORD `.pb/.pb.gz`, CSV, JSON, and custom formats).
- A **standard ReactionRecord dataclass** to hold structured reaction information, conditions, yields, and arbitrary metadata.
- A suite of **filtering and cleaning utilities** to remove noisy entries, canonicalize SMILES, and encode domain heuristics.
- **Reporting tools** that summarize dataset statistics before and after cleaning.
- Optional **PyTorch dataset utilities** for rapid experimentation in reaction-prediction models.

These capabilities allow researchers to focus on modeling and analysis rather than data wrangling, which is often a major barrier to reproducibility in reaction ML projects.

# State of the Field

Reaction-ML workflows in academia and industry rely heavily on data sources such as the USPTO patents dataset, curated proprietary corpora, or structured records from the Open Reaction Database. While several libraries provide downstream modeling components (e.g., reaction prediction architectures, graph neural networks), few open-source tools directly address the *data preparation* stage. Most publications either write custom scripts or use ad hoc pipelines that are rarely reusable across datasets.

chemrxn-cleaner complements existing chemistry toolkits (e.g., RDKit, ord-schema) by focusing explicitly on dataset ingestion, cleaning, and transformation. Its registry-based loader system enables flexible support for new formats, and its filter stack design allows users to encode domain-specific rules—for example, excluding reactions with forbidden elements, filtering by SMILES length, or applying similarity-based criteria. The package therefore helps standardize workflows that are currently duplicated across many research groups.

# Software Description

## Core Features

### 1. Unified data loading
The `load_reactions` function dispatches to format-specific loaders based on a registry. Built-in formats include:

- **USPTO `.rsmi` files**, with optional metadata preservation.
- **ORD `.pb/.pb.gz` files**, extracting identifiers, conditions, yields, and reaction indices.
- **CSV and JSON**, with flexible column mapping and user-provided mappers.
- **Custom formats**, registered with `register_input_format`.

Each loader produces a list of `ReactionRecord` objects.

### 2. Structured reaction representation
`ReactionRecord` stores:

- Parsed SMILES fields (`reaction_smiles`, `reactants`, `reagents`, `products`)
- Optional reaction conditions (temperature, time, solvents, catalysts, additives, pressure, pH, atmosphere, scale)
- Yield information
- Identifiers and provenance
- Arbitrary metadata through `extra_metadata`

Records can be serialized via `to_dict` / `from_dict` and visualized with RDKit.

### 3. Cleaning and filtering
The package includes:

- `clean_reactions` and `clean_and_canonicalize` for parsing and canonicalization.
- Built-in filters such as:
  - `has_product`
  - `all_molecules_valid`
  - `max_smiles_length`
  - `element_filter`
  - `meta_filter`
  - RDKit-based `similarity_filter`

Filters are simple callables returning booleans and may be freely composed.

### 4. Reporting and exporting
The `reporter` module exposes the `CleaningStats`/`FilterStats` structures used by the cleaning pipeline to track how many reactions were kept, dropped, or failed to parse. Cleaned reactions can be exported to JSON or CSV via `export_reaction_records`.

### 5. ML utilities
For rapid experimentation, the package includes:

- `records_to_dataframe` for exploratory analysis.
- `train_valid_test_split` for deterministic dataset partitioning.
- `ForwardReactionDataset`, a minimal PyTorch dataset for forward prediction tasks.

# Example

```python
from chemrxn_cleaner import (
    load_reactions,
    clean_reactions_with_report,
    export_reaction_records,
)

raw = load_reactions("data/sample.rsmi", fmt="uspto", keep_meta=True)

cleaned, stats = clean_reactions_with_report(raw)
print(f"Kept {stats.n_output}/{stats.n_input} reactions after cleaning")

export_reaction_records(cleaned, "cleaned.json", fmt="json")
export_reaction_records(cleaned, "cleaned.csv", fmt="csv")
```
# Acknowledgements
The package builds upon RDKit for cheminformatics operations and ord-schema for parsing ORD records. We thank contributors from the open-source chemistry community whose tools made this project possible.


# References
