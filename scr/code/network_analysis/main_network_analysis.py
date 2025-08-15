#!/usr/bin/env python3
"""
Main Network Analysis Script for TEM Lactamase Proteins

This script orchestrates comprehensive network analysis for TEM beta-lactamase proteins,
comparing results with and without synthetic organisms.

The analysis includes:
- Degree distribution analysis
- Betweenness centrality analysis  
- Graph motif analysis
- Long distance interference analysis

Each analysis is run twice:
1. Excluding synthetic organisms (taxonomy_id = 32630)
2. Including all organisms

Note: All analyses are performed on proteins of length 286 only.

Usage:
    nohup python scr/code/network_analysis/main_network_analysis.py > network_analysis.log 2>&1 &
"""

import logging
import sys
from pathlib import Path
import pandas as pd
import networkx as nx

# =============================================================================
# CONFIGURATION PARAMETERS - MODIFY THESE AS NEEDED
# =============================================================================

# Output directory for all results
OUTPUT_DIR = '/home/nab/Niklas/TEM-lactamase/data/001_results/011_Staan_Stats'

# Standard numbering tool to use for analysis
STANDARD_NUMBERING_TOOL = 'standard_numbering_pairwise_flagged_proteins'

# Analysis options - set to True/False to control which analyses to run
RUN_EXCLUDING_SYNTHETIC = True  # Run analysis excluding synthetic organisms
RUN_INCLUDING_SYNTHETIC = True  # Run analysis including synthetic organisms

# Additional network configurations
RUN_MULTI_MUTATIONS = True      # Include proteins connected by multiple mutations
RUN_ALL_NODES = False          # Include all proteins regardless of sequence length (DISABLED)
RUN_HETEROGENEOUS_NETWORK = True  # Include heterogeneous network with multiple node types

# Add the current directory to Python path for imports
sys.path.append(str(Path(__file__).parent))

from network_core import (
    setup_logging, 
    setup_database_connection, 
    create_output_directory,
    get_tem_family_data,
    get_protein_network_data,
    create_networkx_graph,
    get_heterogeneous_network_data,
    create_heterogeneous_networkx_graph,
    analyze_database_schema,
    save_network_info
)
from network_degree_analysis import analyze_node_degrees
from network_centrality_analysis import analyze_betweenness_centrality
from network_motif_analysis import analyze_graph_motifs
from network_distance_analysis import analyze_long_distance_interference


def run_complete_analysis(eedb, data_tem_ids, output_path, exclude_synthetic=True, allow_multi_mutations=False, include_all_lengths=False, name_of_standard_numbering_tool="standard_numbering_pairwise_flagged_proteins"):
    """Run complete network analysis for a given configuration
    
    Parameters:
    -----------
    eedb : Pyeed
        Database connection object
    data_tem_ids : dict
        TEM family data mapping
    output_path : str
        Directory for saving outputs
    exclude_synthetic : bool
        Whether to exclude synthetic organisms
    allow_multi_mutations : bool
        Whether to include proteins connected by multiple mutations
    include_all_lengths : bool
        Whether to include proteins of all sequence lengths
    name_of_standard_numbering_tool : str
        Standard numbering tool name
        
    Returns:
    --------
    dict : Analysis results summary
    """
    LOGGER = logging.getLogger(__name__)
    
    # Build descriptive analysis type
    organism_type = "excluding synthetic organisms" if exclude_synthetic else "including synthetic organisms"
    mutation_type = "multi-mutations" if allow_multi_mutations else "single mutations"
    length_type = "all lengths" if include_all_lengths else "length 286"
    
    analysis_description = f"{organism_type}, {mutation_type}, {length_type}"
    LOGGER.info(f"Starting analysis: {analysis_description}")
    
    # Build file suffix for outputs
    suffix_parts = []
    suffix_parts.append("excl_synthetic" if exclude_synthetic else "incl_synthetic")
    if allow_multi_mutations:
        suffix_parts.append("multi_mut")
    if include_all_lengths:
        suffix_parts.append("all_lengths")
    
    file_suffix = "_".join(suffix_parts)
    
    # Get protein network data
    LOGGER.info(f"Fetching protein network data...")
    results = get_protein_network_data(eedb, name_of_standard_numbering_tool, exclude_synthetic, allow_multi_mutations, include_all_lengths)
    
    if not results:
        LOGGER.warning(f"No network data found for analysis: {analysis_description}")
        return {}
    
    # Create NetworkX graph
    LOGGER.info(f"Creating NetworkX graph...")
    G = create_networkx_graph(results)
    
    # Save basic network info
    with open(f"{output_path}/network_info_protein_{file_suffix}.txt", 'w') as f:
        f.write(f"NETWORK INFORMATION: PROTEIN NETWORK\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Analysis configuration: {analysis_description}\n")
        f.write(f"Number of nodes: {G.number_of_nodes()}\n")
        f.write(f"Number of edges: {G.number_of_edges()}\n")
        f.write(f"Network density: {nx.density(G):.6f}\n")
        f.write(f"Is connected: {nx.is_connected(G)}\n")
        f.write(f"Number of connected components: {nx.number_connected_components(G)}\n")
        
        if G.number_of_nodes() > 0:
            largest_cc_size = len(max(nx.connected_components(G), key=len))
            f.write(f"Largest connected component size: {largest_cc_size}\n")
            f.write(f"Largest component fraction: {largest_cc_size / G.number_of_nodes():.3f}\n")
    
    # Run all analyses
    analysis_results = {}
    
    LOGGER.info(f"1. Analyzing node degrees...")
    try:
        # For now, use exclude_synthetic but we'll improve this later
        degree_stats = analyze_node_degrees(G, output_path, exclude_synthetic)
        analysis_results['degree_stats'] = degree_stats
    except Exception as e:
        LOGGER.error(f"Degree analysis failed: {e}")
        analysis_results['degree_stats'] = {}
    
    LOGGER.info(f"2. Analyzing betweenness centrality...")
    try:
        centrality_stats = analyze_betweenness_centrality(G, output_path, data_tem_ids, exclude_synthetic)
        analysis_results['centrality_stats'] = centrality_stats
    except Exception as e:
        LOGGER.error(f"Centrality analysis failed: {e}")
        analysis_results['centrality_stats'] = {}
    
    LOGGER.info(f"3. Analyzing graph motifs...")
    try:
        motif_stats = analyze_graph_motifs(G, output_path, exclude_synthetic)
        analysis_results['motif_stats'] = motif_stats
    except Exception as e:
        LOGGER.error(f"Motif analysis failed: {e}")
        analysis_results['motif_stats'] = {}
    
    LOGGER.info(f"4. Analyzing long distance interference...")
    try:
        distance_stats = analyze_long_distance_interference(G, output_path, exclude_synthetic)
        analysis_results['distance_stats'] = distance_stats
    except Exception as e:
        LOGGER.error(f"Distance analysis failed: {e}")
        analysis_results['distance_stats'] = {}
    
    # Mutation length analysis removed as requested
    
    # Store network info
    analysis_results['network_info'] = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'density': G.density() if hasattr(G, 'density') else 0,
        'connected': G.is_connected() if hasattr(G, 'is_connected') else False
    }
    
    LOGGER.info(f"Analysis completed ({organism_type}): {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return analysis_results


def run_heterogeneous_network_analysis(eedb, data_tem_ids, output_path, exclude_synthetic=True):
    """Run complete analysis on heterogeneous network including multiple node types
    
    Parameters:
    -----------
    eedb : Pyeed
        Database connection object
    data_tem_ids : dict
        TEM family data mapping
    output_path : str
        Directory for saving outputs
    exclude_synthetic : bool
        Whether to exclude synthetic organisms
        
    Returns:
    --------
    dict : Analysis results summary
    """
    LOGGER = logging.getLogger(__name__)
    
    # Build descriptive analysis type
    organism_type = "excluding synthetic organisms" if exclude_synthetic else "including synthetic organisms"
    
    analysis_description = f"heterogeneous network, {organism_type}"
    LOGGER.info(f"Starting heterogeneous network analysis: {analysis_description}")
    
    # Build file suffix for outputs
    suffix_parts = []
    suffix_parts.append("excl_synthetic" if exclude_synthetic else "incl_synthetic")
    suffix_parts.append("heterogeneous")
    
    file_suffix = "_".join(suffix_parts)
    
    # Get heterogeneous network data
    LOGGER.info(f"Fetching heterogeneous network data...")
    relationships = get_heterogeneous_network_data(eedb, exclude_synthetic)
    
    if not relationships:
        LOGGER.warning(f"No heterogeneous network data found for analysis: {analysis_description}")
        return {}
    
    # Create NetworkX graph
    LOGGER.info(f"Creating heterogeneous NetworkX graph...")
    G = create_heterogeneous_networkx_graph(relationships)
    
    # Analyze node types
    node_types = {}
    for node, data in G.nodes(data=True):
        node_type = data.get('node_type', 'Unknown')
        node_types[node_type] = node_types.get(node_type, 0) + 1
    
    # Save basic network info
    with open(f"{output_path}/network_info_heterogeneous_{file_suffix}.txt", 'w') as f:
        f.write(f"NETWORK INFORMATION: HETEROGENEOUS NETWORK\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Analysis configuration: {analysis_description}\n")
        f.write(f"Total number of nodes: {G.number_of_nodes()}\n")
        f.write(f"Total number of edges: {G.number_of_edges()}\n")
        f.write(f"Network density: {nx.density(G):.6f}\n")
        f.write(f"Is connected: {nx.is_connected(G)}\n")
        f.write(f"Number of connected components: {nx.number_connected_components(G)}\n")
        
        f.write(f"\nNODE TYPES:\n")
        for node_type, count in node_types.items():
            f.write(f"  {node_type}: {count} nodes\n")
        
        # Analyze edge types
        edge_types = {}
        for u, v, data in G.edges(data=True):
            edge_type = data.get('relationship_type', 'Unknown')
            edge_types[edge_type] = edge_types.get(edge_type, 0) + 1
        
        f.write(f"\nEDGE TYPES:\n")
        for edge_type, count in edge_types.items():
            f.write(f"  {edge_type}: {count} edges\n")
        
        if G.number_of_nodes() > 0:
            largest_cc_size = len(max(nx.connected_components(G), key=len))
            f.write(f"\nLargest connected component size: {largest_cc_size}\n")
            f.write(f"Largest component fraction: {largest_cc_size / G.number_of_nodes():.3f}\n")
    
    # Run analyses (adapted for heterogeneous network)
    analysis_results = {}
    
    LOGGER.info(f"1. Analyzing node degrees (heterogeneous)...")
    try:
        degree_stats = analyze_node_degrees(G, output_path, exclude_synthetic, suffix_override=file_suffix)
        analysis_results['degree_stats'] = degree_stats
    except Exception as e:
        LOGGER.error(f"Degree analysis failed: {e}")
        analysis_results['degree_stats'] = {}
    
    LOGGER.info(f"2. Analyzing betweenness centrality (heterogeneous)...")
    try:
        centrality_stats = analyze_betweenness_centrality(G, output_path, data_tem_ids, exclude_synthetic, suffix_override=file_suffix)
        analysis_results['centrality_stats'] = centrality_stats
    except Exception as e:
        LOGGER.error(f"Centrality analysis failed: {e}")
        analysis_results['centrality_stats'] = {}
    
    LOGGER.info(f"3. Analyzing graph motifs (heterogeneous)...")
    try:
        motif_stats = analyze_graph_motifs(G, output_path, exclude_synthetic, suffix_override=file_suffix)
        analysis_results['motif_stats'] = motif_stats
    except Exception as e:
        LOGGER.error(f"Motif analysis failed: {e}")
        analysis_results['motif_stats'] = {}
    
    LOGGER.info(f"4. Analyzing long distance interference (heterogeneous)...")
    try:
        distance_stats = analyze_long_distance_interference(G, output_path, exclude_synthetic, suffix_override=file_suffix)
        analysis_results['distance_stats'] = distance_stats
    except Exception as e:
        LOGGER.error(f"Distance analysis failed: {e}")
        analysis_results['distance_stats'] = {}
    
    # Store network info
    analysis_results['network_info'] = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'density': nx.density(G),
        'connected': nx.is_connected(G),
        'node_types': node_types,
        'edge_types': edge_types
    }
    
    LOGGER.info(f"Heterogeneous network analysis completed ({organism_type}): {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return analysis_results


def compare_synthetic_vs_natural(results_excl, results_incl, output_path):
    """Compare analysis results between synthetic-excluded and synthetic-included networks
    
    Parameters:
    -----------
    results_excl : dict
        Analysis results excluding synthetic organisms
    results_incl : dict
        Analysis results including synthetic organisms
    output_path : str
        Directory for saving outputs
    """
    LOGGER = logging.getLogger(__name__)
    
    LOGGER.info("Creating comparison between synthetic-excluded and synthetic-included networks...")
    
    # Create comparison summary
    with open(f"{output_path}/synthetic_organism_comparison.txt", 'w') as f:
        f.write("SYNTHETIC ORGANISM IMPACT ANALYSIS\n")
        f.write("=================================\n\n")
        f.write("This analysis compares network properties when synthetic organisms\n")
        f.write("(taxonomy_id = 32630) are included vs excluded from the analysis.\n\n")
        
        # Network size comparison
        f.write("NETWORK SIZE COMPARISON:\n")
        f.write("=======================\n")
        
        excl_nodes = results_excl.get('network_info', {}).get('nodes', 0)
        incl_nodes = results_incl.get('network_info', {}).get('nodes', 0)
        excl_edges = results_excl.get('network_info', {}).get('edges', 0)
        incl_edges = results_incl.get('network_info', {}).get('edges', 0)
        
        f.write(f"Excluding synthetic organisms:\n")
        f.write(f"  Nodes: {excl_nodes}\n")
        f.write(f"  Edges: {excl_edges}\n")
        f.write(f"Including synthetic organisms:\n")
        f.write(f"  Nodes: {incl_nodes}\n")
        f.write(f"  Edges: {incl_edges}\n")
        
        if excl_nodes > 0:
            node_increase = ((incl_nodes - excl_nodes) / excl_nodes) * 100
            edge_increase = ((incl_edges - excl_edges) / excl_edges) * 100 if excl_edges > 0 else 0
            f.write(f"Increase with synthetic organisms:\n")
            f.write(f"  Nodes: +{node_increase:.1f}%\n")
            f.write(f"  Edges: +{edge_increase:.1f}%\n\n")
        
        # Degree statistics comparison
        f.write("DEGREE STATISTICS COMPARISON:\n")
        f.write("============================\n")
        
        excl_degree = results_excl.get('degree_stats', {})
        incl_degree = results_incl.get('degree_stats', {})
        
        if excl_degree and incl_degree:
            f.write(f"Mean degree (excluding synthetic): {excl_degree.get('Mean', 0):.3f}\n")
            f.write(f"Mean degree (including synthetic): {incl_degree.get('Mean', 0):.3f}\n")
            f.write(f"Max degree (excluding synthetic): {excl_degree.get('Max', 0)}\n")
            f.write(f"Max degree (including synthetic): {incl_degree.get('Max', 0)}\n\n")
        
        # Distance statistics comparison
        f.write("DISTANCE STATISTICS COMPARISON:\n")
        f.write("==============================\n")
        
        excl_distance = results_excl.get('distance_stats', {})
        incl_distance = results_incl.get('distance_stats', {})
        
        if excl_distance and incl_distance:
            f.write(f"Diameter (excluding synthetic): {excl_distance.get('diameter', 0)}\n")
            f.write(f"Diameter (including synthetic): {incl_distance.get('diameter', 0)}\n")
            f.write(f"Avg path length (excluding synthetic): {excl_distance.get('avg_path_length', 0):.3f}\n")
            f.write(f"Avg path length (including synthetic): {incl_distance.get('avg_path_length', 0):.3f}\n")
            f.write(f"Small-world sigma (excluding synthetic): {excl_distance.get('small_world_sigma', 0):.3f}\n")
            f.write(f"Small-world sigma (including synthetic): {incl_distance.get('small_world_sigma', 0):.3f}\n\n")
        
        # Motif comparison
        f.write("MOTIF STATISTICS COMPARISON:\n")
        f.write("===========================\n")
        
        excl_motifs = results_excl.get('motif_stats', {})
        incl_motifs = results_incl.get('motif_stats', {})
        
        if excl_motifs and incl_motifs:
            for motif_type in set(list(excl_motifs.keys()) + list(incl_motifs.keys())):
                excl_count = excl_motifs.get(motif_type, 0)
                incl_count = incl_motifs.get(motif_type, 0)
                f.write(f"{motif_type}:\n")
                f.write(f"  Excluding synthetic: {excl_count}\n")
                f.write(f"  Including synthetic: {incl_count}\n")
                if excl_count > 0:
                    increase = ((incl_count - excl_count) / excl_count) * 100
                    f.write(f"  Change: +{increase:.1f}%\n")
                f.write("\n")
        
        # Key findings
        f.write("KEY FINDINGS:\n")
        f.write("============\n")
        
        if incl_nodes > excl_nodes:
            f.write("✓ Synthetic organisms contribute additional proteins to the network\n")
        
        if incl_edges > excl_edges:
            f.write("✓ Synthetic organisms create additional mutation connections\n")
        
        # Check if top centrality nodes changed
        f.write("\nRECOMMENDATION:\n")
        f.write("==============\n")
        
        if incl_nodes > excl_nodes * 1.2:  # More than 20% increase
            f.write("→ Analyze both networks separately as synthetic organisms\n")
            f.write("  significantly affect network structure\n")
        elif incl_nodes > excl_nodes * 1.05:  # More than 5% increase
            f.write("→ Compare top centrality nodes between both networks\n")
            f.write("  to assess biological relevance\n")
        else:
            f.write("→ Synthetic organisms have minimal impact on network structure\n")


def main():
    """Main function to orchestrate complete network analysis"""
    # Setup
    LOGGER = setup_logging()
    LOGGER.info("Starting comprehensive network analysis with synthetic organism comparison")
    
    output_path = create_output_directory(OUTPUT_DIR)
    eedb = setup_database_connection()
    
    # Analyze database schema
    LOGGER.info("Analyzing database schema...")
    schema_info = analyze_database_schema(eedb)
    
    # Get TEM family data
    LOGGER.info("Fetching TEM family data...")
    data_tem_ids = get_tem_family_data(eedb)
    
    # Use the parameters defined at the top of the file
    run_excluding = RUN_EXCLUDING_SYNTHETIC
    run_including = RUN_INCLUDING_SYNTHETIC
    
    results = {}
    
    # Run analysis excluding synthetic organisms
    if run_excluding:
        LOGGER.info("="*60)
        LOGGER.info("ANALYSIS 1: EXCLUDING SYNTHETIC ORGANISMS (SINGLE MUTATIONS, LENGTH 286)")
        LOGGER.info("="*60)
        
        results['excluding_synthetic'] = run_complete_analysis(
            eedb, data_tem_ids, output_path, 
            exclude_synthetic=True, 
            allow_multi_mutations=False,
            include_all_lengths=False,
            name_of_standard_numbering_tool=STANDARD_NUMBERING_TOOL
        )
    
    # Run analysis including synthetic organisms
    if run_including:
        LOGGER.info("="*60)
        LOGGER.info("ANALYSIS 2: INCLUDING SYNTHETIC ORGANISMS (SINGLE MUTATIONS, LENGTH 286)")
        LOGGER.info("="*60)
        
        results['including_synthetic'] = run_complete_analysis(
            eedb, data_tem_ids, output_path, 
            exclude_synthetic=False, 
            allow_multi_mutations=False,
            include_all_lengths=False,
            name_of_standard_numbering_tool=STANDARD_NUMBERING_TOOL
        )
    
    # Run analysis with multi-mutations (excluding synthetic organisms)
    if RUN_MULTI_MUTATIONS:
        LOGGER.info("="*60)
        LOGGER.info("ANALYSIS 3: MULTI-MUTATIONS NETWORK (EXCLUDING SYNTHETIC, LENGTH 286)")
        LOGGER.info("="*60)
        
        results['multi_mutations'] = run_complete_analysis(
            eedb, data_tem_ids, output_path, 
            exclude_synthetic=True, 
            allow_multi_mutations=True,
            include_all_lengths=False,
            name_of_standard_numbering_tool=STANDARD_NUMBERING_TOOL
        )
    
    # All nodes network analysis disabled as requested
    if RUN_ALL_NODES:
        pass  # All nodes network analysis disabled as requested
    
    # Run heterogeneous network analysis (excluding synthetic organisms)
    if RUN_HETEROGENEOUS_NETWORK:
        LOGGER.info("="*60)
        LOGGER.info("ANALYSIS 4: HETEROGENEOUS NETWORK (EXCLUDING SYNTHETIC)")
        LOGGER.info("="*60)
        
        results['heterogeneous_network'] = run_heterogeneous_network_analysis(
            eedb, data_tem_ids, output_path, 
            exclude_synthetic=True
        )
    
    # Create comparison if both analyses were run
    if run_excluding and run_including:
        LOGGER.info("="*60)
        LOGGER.info("COMPARISON: SYNTHETIC VS NATURAL ORGANISMS")
        LOGGER.info("="*60)
        
        compare_synthetic_vs_natural(
            results['excluding_synthetic'], 
            results['including_synthetic'], 
            output_path
        )
    
    # Create comprehensive summary
    LOGGER.info("Creating comprehensive summary...")
    
    with open(f"{output_path}/comprehensive_analysis_summary.txt", 'w') as f:
        f.write("COMPREHENSIVE TEM LACTAMASE NETWORK ANALYSIS SUMMARY\n")
        f.write("===================================================\n\n")
        f.write(f"Analysis completed on: {pd.Timestamp.now()}\n")
        f.write(f"Output directory: {output_path}\n")
        f.write(f"Standard numbering tool: {STANDARD_NUMBERING_TOOL}\n\n")
        
        f.write("DATABASE SCHEMA OVERVIEW:\n")
        f.write(f"Available node types: {len(schema_info['node_labels'])}\n")
        for label in schema_info['node_labels']:
            count = schema_info['node_counts'].get(label, 0)
            f.write(f"  - {label}: {count} nodes\n")
        f.write(f"Available relationship types: {len(schema_info['relationship_types'])}\n")
        for rel_type in schema_info['relationship_types']:
            f.write(f"  - {rel_type}\n")
        f.write("\n")
        
        # Analysis results summary
        if 'excluding_synthetic' in results:
            excl_net = results['excluding_synthetic'].get('network_info', {})
            f.write("ANALYSIS 1 - EXCLUDING SYNTHETIC ORGANISMS (single mutations, length 286):\n")
            f.write(f"  Nodes: {excl_net.get('nodes', 0)}\n")
            f.write(f"  Edges: {excl_net.get('edges', 0)}\n")
            f.write(f"  Density: {excl_net.get('density', 0):.6f}\n\n")
        
        if 'including_synthetic' in results:
            incl_net = results['including_synthetic'].get('network_info', {})
            f.write("ANALYSIS 2 - INCLUDING SYNTHETIC ORGANISMS (single mutations, length 286):\n")
            f.write(f"  Nodes: {incl_net.get('nodes', 0)}\n")
            f.write(f"  Edges: {incl_net.get('edges', 0)}\n")
            f.write(f"  Density: {incl_net.get('density', 0):.6f}\n\n")
        
        if 'multi_mutations' in results:
            multi_net = results['multi_mutations'].get('network_info', {})
            f.write("ANALYSIS 3 - MULTI-MUTATIONS NETWORK (excluding synthetic, length 286):\n")
            f.write(f"  Nodes: {multi_net.get('nodes', 0)}\n")
            f.write(f"  Edges: {multi_net.get('edges', 0)}\n")
            f.write(f"  Density: {multi_net.get('density', 0):.6f}\n\n")
        
        # All nodes analysis removed as requested
        
        if 'heterogeneous_network' in results:
            hetero_net = results['heterogeneous_network'].get('network_info', {})
            f.write("ANALYSIS 4 - HETEROGENEOUS NETWORK (excluding synthetic):\n")
            f.write(f"  Nodes: {hetero_net.get('nodes', 0)}\n")
            f.write(f"  Edges: {hetero_net.get('edges', 0)}\n")
            f.write(f"  Density: {hetero_net.get('density', 0):.6f}\n")
            node_types = hetero_net.get('node_types', {})
            if node_types:
                f.write(f"  Node types:\n")
                for node_type, count in node_types.items():
                    f.write(f"    {node_type}: {count}\n")
            f.write("\n")
        
        f.write("FILES GENERATED:\n")
        f.write("===============\n")
        if run_excluding:
            f.write("ANALYSIS 1 - EXCLUDING SYNTHETIC (single mutations, length 286):\n")
            f.write("  - degree_analysis_excl_synthetic.png\n")
            f.write("  - betweenness_centrality_analysis_excl_synthetic.png\n")
            f.write("  - motif_analysis_excl_synthetic.png\n")
            f.write("  - long_distance_interference_analysis_excl_synthetic.png\n")
            f.write("  - Various statistics files with _excl_synthetic suffix\n\n")
        
        if run_including:
            f.write("ANALYSIS 2 - INCLUDING SYNTHETIC (single mutations, length 286):\n")
            f.write("  - degree_analysis_incl_synthetic.png\n")
            f.write("  - betweenness_centrality_analysis_incl_synthetic.png\n")
            f.write("  - motif_analysis_incl_synthetic.png\n")
            f.write("  - long_distance_interference_analysis_incl_synthetic.png\n")
            f.write("  - Various statistics files with _incl_synthetic suffix\n\n")
        
        if RUN_MULTI_MUTATIONS:
            f.write("ANALYSIS 3 - MULTI-MUTATIONS (excluding synthetic, length 286):\n")
            f.write("  - degree_analysis_excl_synthetic_multi_mut.png\n")
            f.write("  - betweenness_centrality_analysis_excl_synthetic_multi_mut.png\n")
            f.write("  - motif_analysis_excl_synthetic_multi_mut.png\n")
            f.write("  - long_distance_interference_analysis_excl_synthetic_multi_mut.png\n")
            f.write("  - Various statistics files with _excl_synthetic_multi_mut suffix\n\n")
        
        # All nodes analysis removed as requested
        
        if RUN_HETEROGENEOUS_NETWORK:
            f.write("ANALYSIS 4 - HETEROGENEOUS NETWORK (excluding synthetic):\n")
            f.write("  - degree_analysis_excl_synthetic_heterogeneous.png\n")
            f.write("  - betweenness_centrality_analysis_excl_synthetic_heterogeneous.png\n")
            f.write("  - motif_analysis_excl_synthetic_heterogeneous.png\n")
            f.write("  - long_distance_interference_analysis_excl_synthetic_heterogeneous.png\n")
            f.write("  - Various statistics files with _excl_synthetic_heterogeneous suffix\n\n")
        
        if run_excluding and run_including:
            f.write("COMPARISON FILES:\n")
            f.write("  - synthetic_organism_comparison.txt\n")
    
    print(f"\nComprehensive network analysis completed!")
    print(f"Results saved to: {output_path}")
    print(f"\nNetwork Analysis Summary:")
    
    if 'excluding_synthetic' in results:
        excl_net = results['excluding_synthetic']['network_info']
        print(f"1. Excluding synthetic (single mutations, length 286): {excl_net['nodes']} nodes, {excl_net['edges']} edges")
    
    if 'including_synthetic' in results:
        incl_net = results['including_synthetic']['network_info']
        print(f"2. Including synthetic (single mutations, length 286): {incl_net['nodes']} nodes, {incl_net['edges']} edges")
    
    if 'multi_mutations' in results:
        multi_net = results['multi_mutations']['network_info']
        print(f"3. Multi-mutations (excluding synthetic, length 286): {multi_net['nodes']} nodes, {multi_net['edges']} edges")
    
    # All nodes analysis removed as requested
    
    if 'heterogeneous_network' in results:
        hetero_net = results['heterogeneous_network']['network_info']
        node_types_summary = ", ".join([f"{t}: {c}" for t, c in hetero_net.get('node_types', {}).items()])
        print(f"4. Heterogeneous network (excluding synthetic): {hetero_net['nodes']} nodes, {hetero_net['edges']} edges ({node_types_summary})")
    
    if run_excluding and run_including:
        print("\nComparison analysis completed - check synthetic_organism_comparison.txt for key differences")
    
    total_analyses = len([k for k in results.keys() if k in ['excluding_synthetic', 'including_synthetic', 'multi_mutations', 'heterogeneous_network']])
    print(f"\nTotal of {total_analyses} network configurations analyzed.")
    
    LOGGER.info(f"Analysis completed successfully. Results saved to: {output_path}")


if __name__ == "__main__":
    main()
