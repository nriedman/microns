"""
Microns - Synaptic Organization Analysis Package

This package contains modules for analyzing synaptic organization in dendrites
using minimum spanning trees from the MICrONS Cortical MM^3 dataset.

Modules:
    mst_generation: Core functions for generating minimum spanning trees
    visualization: Functions for plotting and visualizing MSTs
    examples: Usage examples and workflows
"""

from .mst_generation import generate_msts, connect_disjoint_branches
from .visualization import plot_mst_3d, traverse_branch

__all__ = ['generate_msts', 'connect_disjoint_branches', 'plot_mst_3d', 'traverse_branch']