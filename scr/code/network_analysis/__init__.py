"""
TEM Lactamase Network Analysis Package

This package provides modular network analysis tools for TEM beta-lactamase proteins.
It supports analysis both with and without synthetic organisms to compare their impact
on network structure and properties.

Modules:
--------
- network_core: Core functions for database connection, data fetching, and graph creation
- network_degree_analysis: Node degree distribution analysis
- network_centrality_analysis: Betweenness centrality analysis
- network_motif_analysis: Graph motif and structural pattern analysis
- network_distance_analysis: Path length and small-world analysis
- network_phenotype_analysis: Phenotype-based mutation analysis
- main_network_analysis: Main orchestrator script

Usage:
------
Run complete analysis:
    python main_network_analysis.py

Run with custom options:
    python main_network_analysis.py --output-dir /path/to/output --exclude-synthetic-only

Import individual modules:
    from network_analysis.network_degree_analysis import analyze_node_degrees
    from network_analysis.network_core import setup_database_connection
"""

__version__ = "1.0.0"
__author__ = "Network Analysis Team"
