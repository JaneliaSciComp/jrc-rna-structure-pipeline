## Overview

This pipeline automates the extraction and processing of RNA structures from PDB, including:

- Fetching RNA structures within a specified date range
- Extracting metadata and ligand information (SMILES)
- Grouping and clustering sequences to remove redundancy
- Calculating structural metrics (structuredness)
- Filtering based on quality criteria
- Splitting datasets into training and testing sets
- Converting to Kaggle competition format

## Prerequisites

- [Pixi](https://pixi.sh/) package manager
- Linux environment (tested on Linux-64)

## Installation

### 1. Install Pixi

If you don't have Pixi installed, follow the instructions at [https://pixi.sh/](https://pixi.sh/):

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### 2. Set up the environment

Clone this repository and navigate to the project directory:

```bash
cd jrc-rna-structure-pipeline
```

Install dependencies using Pixi:

```bash
pixi install
```

This will create a conda environment with all required dependencies including:

- Python packages: pandas, biopython, duckdb
- Bioinformatics tools: mmseqs2

### 3. Activate the environment

```bash
pixi shell
```

## Usage

### Running the Kaggle Workflow

The easiest way to run the pipeline is using one of the workflow scripts used to generate. Here's an example using the workflow used to compile data for the
[Stanford RNA 3D Folding Part 2](https://www.kaggle.com/competitions/stanford-rna-3d-folding-2) competition. Note: workflow script is setup to be run in the output directory.

```bash
cd /path/to/working/directory
bash /path/to/jrc-rna-structure-pipeline/workflows/kaggle2026_1978-01-01_2025-12-17.sh
```

This workflow:

- Fetches RNA structures from PDB (1978-01-01 to 2025-12-17)
- Processes and filters the data
- Creates train/test split at 2025-05-29
- Generates Kaggle-formatted output in `kaggle_jrc_v330/`

### Workflow Options

The workflow script supports various options for controlling execution. The steps can be setup

```bash
# Show help and available options
bash workflows/kaggle2026_1978-01-01_2025-12-17.sh --help

# Dry run to see which steps will execute
bash workflows/kaggle2026_1978-01-01_2025-12-17.sh --dry-run

# Skip specific steps
bash workflows/kaggle2026_1978-01-01_2025-12-17.sh --skip-fetch --skip-cluster

# Run only specific groups
bash workflows/kaggle2026_1978-01-01_2025-12-17.sh --only-kaggle
```

### Pipeline Steps

The complete pipeline consists of these major steps:

1. **Fetch Data**: Download RNA structures from PDB within date range
2. **Extract Metadata**: Parse structure files and extract metadata
3. **Extract Bioassembly**: Process biological assemblies
4. **Extract SMILES**: Extract ligand information
5. **Group Sequences**: Remove 100% identical sequences
6. **Cluster with MMseqs2**: Cluster at various identity thresholds (30%, 50%, 70%, 90%)
7. **Calculate Structuredness**: Compute structural metrics for each chain
8. **Add Structuredness**: Merge structural metrics into metadata
9. **Filter Metadata**: Apply quality filters
10. **Select Representatives**: Choose representative structures from clusters
11. **Extract Coordinates**: Extract 3D coordinates for selected structures
12. **Split Train/Test**: Temporal split based on release date
13. **Convert to Kaggle Format**: Generate final output files

### Output Files

The pipeline generates several output files:

- `jrc_v330_metadata.csv`: Metadata with structuredness metrics
- `jrc_v330_metadata_filtered.csv`: Filtered high-quality structures
- `train_jrc_v330_bioassembly.csv`: Training set metadata
- `test_jrc_v330_bioassembly.csv`: Test set metadata
- `kaggle_jrc_v330/`: Kaggle competition format files
  - `train_sequences.csv`: Training sequences
  - `test_sequences.csv`: Test sequences
  - `validation_sequences.csv`: Validation sequences
  - `train_labels.csv`: Training coordinates
  - `test_labels.csv`: Test coordinates
  - `validation_labels.csv`: Validation coordinates

## Workflow backwards compatibility

Currently we do not maintain backwards compatibility between different versions. Ie. the workflow setup for the current version of the repository may not execute correctly
with newer version of the repository. We use versioned releases to track the compatible versions.
