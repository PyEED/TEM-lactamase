import logging
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from collections import Counter


def analyze_node_degrees(G, output_path, exclude_synthetic=True, suffix_override=None):
    """Perform node degree analysis and create visualizations
    
    Note: Edges represent single mutations only (proteins that differ
    by exactly one amino acid position). Degree = number of proteins
    reachable by single mutations.
    
    Parameters:
    -----------
    G : networkx.Graph
        The protein network graph
    output_path : str
        Directory path for saving outputs
    exclude_synthetic : bool
        Whether synthetic organisms were excluded from the analysis
    suffix_override : str, optional
        Override the default file suffix (for heterogeneous networks)
    """
    LOGGER = logging.getLogger(__name__)
    
    # Calculate degree sequence and distribution
    degree_sequence = [d for _, d in G.degree()]
    degree_counts = Counter(degree_sequence)
    degrees = sorted(degree_counts.keys())
    counts = [degree_counts[d] for d in degrees]

    # Log-binning for better visualization
    max_degree = max(degrees) if degrees else 0
    if max_degree > 1:
        bin_edges = np.unique(2 ** np.arange(0, int(np.log2(max_degree))+1))
    else:
        bin_edges = np.array([0, 1])
    
    binned_degrees = []
    binned_counts = []
    for left, right in zip(bin_edges[:-1], bin_edges[1:]):
        mask = [left <= k < right for k in degree_counts.keys()]
        count = sum([c for k, c, m in zip(degree_counts.keys(),
                                         degree_counts.values(),
                                         mask) if m])
        if count > 0:
            binned_degrees.append(left)
            binned_counts.append(count)

    binned_degrees = np.array(binned_degrees, dtype=float)
    binned_counts = np.array(binned_counts, dtype=float)

    # Calculate statistics first for use in titles
    if degree_sequence:
        mean_degree = np.mean(degree_sequence)
        std_degree = np.std(degree_sequence)
        max_degree = max(degree_sequence)
        min_degree = min(degree_sequence)
        stats_text = f"(μ={mean_degree:.2f}, σ={std_degree:.2f}, min={min_degree}, max={max_degree})"
        degree_stats = {
            'Mean': mean_degree,
            'Median': np.median(degree_sequence),
            'Std': std_degree,
            'Max': max_degree,
            'Min': min_degree
        }
    else:
        stats_text = "(no data)"
        degree_stats = {'Mean': 0, 'Median': 0, 'Std': 0, 'Max': 0, 'Min': 0}

    # Create plots
    if suffix_override:
        organism_suffix = f"_{suffix_override}"
        organism_title = "Heterogeneous Network" if "heterogeneous" in suffix_override else ("Excluding Synthetic Organisms" if exclude_synthetic else "Including Synthetic Organisms")
    else:
        organism_suffix = "_excl_synthetic" if exclude_synthetic else "_incl_synthetic"
        organism_title = "Excluding Synthetic Organisms" if exclude_synthetic else "Including Synthetic Organisms"
    
    plt.figure(figsize=(18, 6))
    
    # Subplot 1: Regular degree distribution
    plt.subplot(1, 3, 1)
    plt.bar(degrees, counts, alpha=0.7)
    plt.title(f"Protein Network: Degree Distribution\n(Single Mutation Edges) - {organism_title}\n{stats_text}")
    plt.xlabel("Degree")
    plt.ylabel("Frequency")
    plt.grid(True, alpha=0.3)

    # Subplot 2: Log-log degree distribution with binning
    plt.subplot(1, 3, 2)
    if degrees and max(degrees) > 0:
        plt.loglog(degrees, counts, "bo-", linewidth=1, markersize=4, alpha=0.6, label='Raw data')
        if len(binned_degrees) > 0:
            plt.loglog(binned_degrees, binned_counts, "rs-", linewidth=2, markersize=8, label='Log-binned')
        plt.grid(True, which="both", ls="--", alpha=0.4)
        plt.xlabel("Degree (log)")
        plt.ylabel("Frequency (log)")
        plt.legend()
    plt.title(f"Protein Network: Degree Distribution (Log-Log Scale)\n(Single Mutation Edges) - {organism_title}\n{stats_text}")

    # Subplot 3: Cumulative degree distribution
    plt.subplot(1, 3, 3)
    if degree_sequence:
        sorted_degrees = sorted(degree_sequence, reverse=True)
        cumulative = np.arange(1, len(sorted_degrees) + 1) / len(sorted_degrees)
        if max(sorted_degrees) > 0:
            plt.loglog(sorted_degrees, cumulative, "go-", linewidth=2, markersize=4)
        plt.grid(True, which="both", ls="--", alpha=0.4)
        plt.xlabel("Degree (log)")
        plt.ylabel("P(X ≥ k)")
    plt.title(f"Protein Network: Cumulative Degree Distribution\n(Single Mutation Edges) - {organism_title}\n{stats_text}")

    # Statistics subplot removed - statistics now included in titles

    plt.tight_layout()
    plt.savefig(f"{output_path}/degree_analysis{organism_suffix}.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Save degree statistics to file
    with open(f"{output_path}/degree_statistics{organism_suffix}.txt", 'w') as f:
        f.write("PROTEIN NETWORK: NODE DEGREE ANALYSIS\n")
        f.write("====================================\n\n")
        f.write(f"ORGANISM INCLUSION: {organism_title}\n\n")
        f.write("EDGE CRITERIA:\n")
        f.write("Edges represent proteins connected by single mutations only\n")
        f.write("(proteins that differ by exactly one amino acid position)\n\n")
        f.write("INTERPRETATION:\n")
        f.write("Node degree = number of proteins reachable by single mutations\n")
        f.write("High degree = protein with many single-step mutational neighbors\n\n")
        f.write(f"Number of nodes: {G.number_of_nodes()}\n")
        f.write(f"Number of edges: {G.number_of_edges()}\n")
        
        if degree_sequence:
            f.write(f"Mean degree: {mean_degree:.3f}\n")
            f.write(f"Median degree: {degree_stats['Median']:.3f}\n")
            f.write(f"Standard deviation: {std_degree:.3f}\n")
            f.write(f"Maximum degree: {max_degree}\n")
            f.write(f"Minimum degree: {min_degree}\n")
        else:
            f.write("No degree data available\n")
            
        f.write(f"\nLOG BINNING INFORMATION:\n")
        f.write(f"Number of log bins: {len(binned_degrees)}\n")
        f.write(f"Bin edges: {bin_edges.tolist() if len(bin_edges) > 0 else 'None'}\n")
        f.write(f"Binned degrees: {binned_degrees.tolist() if len(binned_degrees) > 0 else 'None'}\n")
        f.write(f"Binned counts: {binned_counts.tolist() if len(binned_counts) > 0 else 'None'}\n")

    final_mean_degree = mean_degree if degree_sequence else 0
    LOGGER.info(f"Degree analysis completed ({'excl' if exclude_synthetic else 'incl'} synthetic). Mean degree: {final_mean_degree:.3f}")
    return degree_stats if degree_sequence else {}
