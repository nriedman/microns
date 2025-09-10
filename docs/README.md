# MICrONS: Synaptic Organization Analysis

A Summer 2023 research project analyzing the synaptic organization of dendrites in the MICrONS Cortical MM³ dataset using minimum spanning trees (MSTs).

## Overview

This project investigates dendritic organization patterns by generating and analyzing minimum spanning trees from synaptic connectivity data. The approach reveals structural patterns in how synapses are organized along dendrites in cortical neurons.

## Repository Structure

```
microns/
├── src/                    # Core analysis modules
├── notebooks/              # Analysis notebooks (organized by purpose)
├── data/                   # Datasets and mapping files
├── environments/           # Conda environment configurations
└── docs/                   # Documentation
```

## Quick Start

### 1. Environment Setup

```bash
# Clone the repository
git clone <repository-url>
cd microns

# Create conda environment
conda env create -f environments/environment.yml
conda activate microns

# Or use the full environment snapshot
conda env create -f environments/microns_environment.yml
```

### 2. Basic Usage

```python
import pandas as pd
from src.mst_generation import generate_msts
from src.visualization import plot_mst_3d

# Load your data
synapse_df = pd.read_csv('data/sample_datasets/your_synapses.csv', index_col=0)
synapse_df.set_index('synapse_id', inplace=True)

cells_df = pd.read_csv('data/cells_no_repeats.csv', index_col=0)
cells_df.set_index('pt_root_id', inplace=True)

# Generate minimum spanning trees
msts, failed_cells = generate_msts(
    synapses_df=synapse_df,
    cells_df=cells_df,
    k=6,        # nearest neighbors for synapses
    soma_k=10   # nearest neighbors for soma
)

# Visualize the first MST
if msts:
    plot_mst_3d(msts[0])
```

## Core Components

### MST Generation (`src/mst_generation.py`)

The main function `generate_msts()` creates minimum spanning trees from synaptic data:

**Parameters:**
- `synapses_df`: DataFrame of synapses (indexed by synapse_id)
- `cells_df`: DataFrame of cell information (indexed by pt_root_id)  
- `k`: Number of nearest neighbors for regular synapses (default: 6)
- `soma_k`: Number of nearest neighbors for soma connections (default: 10)

**Returns:**
- List of NetworkX Graph objects (MSTs)
- List of cell IDs with insufficient synapses

**Required DataFrame columns:**

*Synapse DataFrame:*
- `ctr_pt_x`, `ctr_pt_y`, `ctr_pt_z`: Synapse coordinates
- `post_pt_root_id`: Postsynaptic cell ID
- `pre_pt_root_id`: Presynaptic cell ID
- `pre_cell_type`: Presynaptic cell type

*Cells DataFrame:*
- `pt_root_id`: Cell ID
- `pt_x`, `pt_y`, `pt_z`: Cell body coordinates
- `cell_type`: Cell type classification

### Visualization (`src/visualization.py`)

Functions for plotting and analyzing MSTs:

- `plot_mst_3d(G)`: Creates 3D visualization with colored sequence paths
- `traverse_branch()`: Extracts linear paths from MST graphs
- Options for burst sequence filtering and cell ID conversion

### Data Organization

- **`data/cells_no_repeats.csv`**: Complete cell information for the MM³ dataset
- **`data/mappings/`**: ID conversion dictionaries for readable sequence representation
- **`data/sample_datasets/`**: Example synapse datasets (numbers indicate input percentile)

## Analysis Notebooks

### Core Analysis (`notebooks/analysis/`)
- `final_analysis.ipynb`: Main results and conclusions
- `minimum_spanning_tree.ipynb`: Core MST methodology
- `minimum_spanning_tree_clean.ipynb`: Cleaned analysis pipeline

### Sequence Analysis (`notebooks/sequences/`)
- `burst_sequences.ipynb`: Synaptic burst pattern analysis
- `extract_sequences_den.ipynb`: Sequence extraction methods

### Evaluation (`notebooks/evaluation/`)
- `test_mst_results.ipynb`: Method validation
- `mst_corner_*.ipynb`: Edge case analysis
- `query_corners.ipynb`: Data exploration

## Advanced Features

### Sequence Analysis

The toolkit includes functionality to extract and analyze synaptic sequences:

```python
from src.visualization import traverse_branch

# Extract sequences from an MST
sequences = []
for start_node in leaf_nodes:
    path = traverse_branch(
        mst_graph, 
        start_node, 
        soma_node, 
        ignore_burst=True,    # Filter burst sequences
        return_seq=True       # Return cell IDs vs synapse IDs
    )
    sequences.append(path)
```

### Cell ID Mapping

Convert between cell IDs and readable character representations:

```python
import pickle

# Load mapping dictionaries
with open('data/mappings/pt_root_id_to_char.pkl', 'rb') as f:
    id_to_char = pickle.load(f)

with open('data/mappings/char_to_codon.pkl', 'rb') as f:
    char_to_codon = pickle.load(f)

# Convert sequences to readable format
readable_sequence = [char_to_codon[id_to_char[cell_id]] for cell_id in sequence]
```

## Key Features

- **Flexible MST Generation**: Configurable nearest neighbor parameters
- **3D Visualization**: Interactive plotting of dendritic organization
- **Sequence Analysis**: Extract and analyze synaptic input patterns
- **Burst Detection**: Optional filtering of repetitive synaptic events
- **ID Mapping**: Convert between numerical IDs and readable formats
- **Modular Design**: Organized codebase for easy extension

## Requirements

- Python 3.8.10
- Core: pandas, numpy, matplotlib, networkx, scikit-learn
- Visualization: jupyterlab, ipympl
- Neuroscience: allensdk, caveclient, meshparty
- See `environments/` for complete dependency lists

## Citation

This work analyzes data from the MICrONS Cortical MM³ dataset. Please cite appropriate MICrONS publications when using this code.

## Questions?

For questions about usage or implementation, please refer to the inline documentation in the source code or the example notebooks.