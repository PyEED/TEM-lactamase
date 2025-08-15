import logging
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict


def analyze_betweenness_centrality(G, output_path, data_tem_ids, exclude_synthetic=True, suffix_override=None):
    """Calculate and analyze betweenness centrality
    
    Note: Edges in this protein network represent single mutations only
    (proteins that differ by exactly one amino acid position).
    Betweenness centrality measures how important each protein is as an
    intermediary in single-step mutation pathways.
    
    Parameters:
    -----------
    G : networkx.Graph
        The protein network graph
    output_path : str
        Directory path for saving outputs
    data_tem_ids : dict
        Mapping of TEM names to protein IDs
    exclude_synthetic : bool
        Whether synthetic organisms were excluded from the analysis
    suffix_override : str, optional
        Override the default file suffix (for heterogeneous networks)
    """
    LOGGER = logging.getLogger(__name__)
    
    # Calculate betweenness centrality
    betweenness_centrality = nx.betweenness_centrality(G)
    
    # Get top nodes by betweenness centrality
    top_nodes = sorted(betweenness_centrality.items(), key=lambda x: x[1], reverse=True)[:20]
    top_nodes_dict = {node: value for node, value in top_nodes}

    # Map node IDs to TEM names where possible
    id_to_name = {}
    for name, ids in data_tem_ids.items():
        for node_id in ids:
            id_to_name[node_id] = name

    labeled_nodes = {}
    for node, value in top_nodes_dict.items():
        if node in id_to_name:
            labeled_nodes[id_to_name[node]] = value
        else:
            labeled_nodes[node[:10] + "..."] = value  # Truncate long IDs

    # Create plots
    if suffix_override:
        organism_suffix = f"_{suffix_override}"
        organism_title = "Heterogeneous Network" if "heterogeneous" in suffix_override else ("Excluding Synthetic Organisms" if exclude_synthetic else "Including Synthetic Organisms")
    else:
        organism_suffix = "_excl_synthetic" if exclude_synthetic else "_incl_synthetic"
        organism_title = "Excluding Synthetic Organisms" if exclude_synthetic else "Including Synthetic Organisms"
    
    # Calculate statistics first for use in titles
    centrality_values = list(betweenness_centrality.values())
    if centrality_values:
        mean_centrality = np.mean(centrality_values)
        std_centrality = np.std(centrality_values)
        max_centrality = max(centrality_values)
        stats_text = f"(μ={mean_centrality:.4f}, σ={std_centrality:.4f}, max={max_centrality:.4f})"
    else:
        stats_text = "(no data)"
    
    plt.figure(figsize=(18, 6))
    
    # Subplot 1: Top betweenness centrality nodes
    plt.subplot(1, 3, 1)
    if labeled_nodes:
        plt.bar(range(len(labeled_nodes)), labeled_nodes.values())
        plt.xticks(range(len(labeled_nodes)), labeled_nodes.keys(), rotation=45, ha='right')
    plt.title(f"Protein Network: Top 20 Nodes by Betweenness Centrality\n(Single Mutation Edges) - {organism_title}\n{stats_text}")
    plt.xlabel("Node")
    plt.ylabel("Betweenness Centrality")

    # Subplot 2: Betweenness centrality distribution
    plt.subplot(1, 3, 2)
    if centrality_values and max(centrality_values) > 0:
        plt.hist(centrality_values, bins=50, alpha=0.7, edgecolor='black')
        plt.yscale('log')
    plt.title(f"Protein Network: Betweenness Centrality Distribution\n(Single Mutation Edges) - {organism_title}\n{stats_text}")
    plt.xlabel("Betweenness Centrality")
    plt.ylabel("Frequency")

    # Subplot 3: Centrality vs Degree correlation
    plt.subplot(1, 3, 3)
    degrees = [G.degree(node) for node in G.nodes()]
    centralities = [betweenness_centrality[node] for node in G.nodes()]
    if degrees and centralities:
        plt.scatter(degrees, centralities, alpha=0.6)
        
        # Calculate correlation
        if len(degrees) > 1 and np.std(degrees) > 0 and np.std(centralities) > 0:
            correlation = np.corrcoef(degrees, centralities)[0, 1]
        else:
            correlation = 0
        plt.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
                 transform=plt.gca().transAxes, fontsize=12, 
                 bbox=dict(boxstyle="round", facecolor='wheat', alpha=0.8))
    else:
        correlation = 0
    plt.xlabel("Node Degree")
    plt.ylabel("Betweenness Centrality")
    plt.title(f"Protein Network: Degree vs Betweenness Centrality\n(Single Mutation Edges) - {organism_title}\n{stats_text}")

    # Statistics subplot removed - statistics now included in titles

    plt.tight_layout()
    plt.savefig(f"{output_path}/betweenness_centrality_analysis{organism_suffix}.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Save centrality statistics
    with open(f"{output_path}/betweenness_centrality_statistics{organism_suffix}.txt", 'w') as f:
        f.write("PROTEIN NETWORK: BETWEENNESS CENTRALITY ANALYSIS\n")
        f.write("===============================================\n\n")
        f.write(f"ORGANISM INCLUSION: {organism_title}\n\n")
        f.write("EDGE CRITERIA:\n")
        f.write("Edges represent proteins connected by single mutations only\n")
        f.write("(proteins that differ by exactly one amino acid position)\n\n")
        f.write("INTERPRETATION:\n")
        f.write("Betweenness centrality measures how important each protein is\n")
        f.write("as an intermediary in single-step mutation pathways.\n")
        f.write("High centrality = key stepping stone in evolutionary paths.\n\n")
        
        if centrality_values:
            f.write(f"Mean betweenness centrality: {mean_centrality:.6f}\n")
            f.write(f"Median betweenness centrality: {np.median(centrality_values):.6f}\n")
            f.write(f"Standard deviation: {std_centrality:.6f}\n")
            f.write(f"Maximum centrality: {max_centrality:.6f}\n")
        else:
            f.write("No centrality data available\n")
            
        f.write(f"Correlation with degree: {correlation:.3f}\n\n")
        
        f.write("TOP 20 NODES BY BETWEENNESS CENTRALITY:\n")
        for i, (node, centrality) in enumerate(top_nodes):
            name = id_to_name.get(node, "Unknown")
            f.write(f"{i+1:2d}. {node} → {name}: {centrality:.6f}\n")

    final_max_centrality = max_centrality if centrality_values else 0
    LOGGER.info(f"Betweenness centrality analysis completed ({'excl' if exclude_synthetic else 'incl'} synthetic). Max centrality: {final_max_centrality:.6f}")
    return betweenness_centrality
