import logging
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from collections import Counter


def analyze_long_distance_interference(G, output_path, exclude_synthetic=True):
    """Investigate long distance interference in the graph
    
    Note: Analyzes path lengths in single-mutation protein network.
    Path length = minimum number of single mutations required to connect two proteins.
    This reveals the mutational distance and evolutionary relationships.
    
    Parameters:
    -----------
    G : networkx.Graph
        The protein network graph
    output_path : str
        Directory path for saving outputs
    exclude_synthetic : bool
        Whether synthetic organisms were excluded from the analysis
    """
    LOGGER = logging.getLogger(__name__)
    
    # Calculate shortest path lengths
    if nx.is_connected(G):
        # For connected graph
        shortest_paths = dict(nx.all_pairs_shortest_path_length(G))
        path_lengths = []
        for source in shortest_paths:
            for target, length in shortest_paths[source].items():
                if source != target:
                    path_lengths.append(length)
        
        diameter = max(path_lengths) if path_lengths else 0
        avg_path_length = np.mean(path_lengths) if path_lengths else 0
        
    else:
        # For disconnected graph, analyze largest component
        if G.number_of_nodes() > 0:
            largest_cc = max(nx.connected_components(G), key=len)
            G_largest = G.subgraph(largest_cc)
            
            shortest_paths = dict(nx.all_pairs_shortest_path_length(G_largest))
            path_lengths = []
            for source in shortest_paths:
                for target, length in shortest_paths[source].items():
                    if source != target:
                        path_lengths.append(length)
            
            diameter = max(path_lengths) if path_lengths else 0
            avg_path_length = np.mean(path_lengths) if path_lengths else 0
        else:
            diameter = 0
            avg_path_length = 0
            path_lengths = []
    
    # Calculate efficiency and resistance distance approximation
    efficiency = 1 / avg_path_length if avg_path_length > 0 else 0
    
    # Analyze path length distribution
    path_length_counts = Counter(path_lengths) if path_lengths else Counter()
    
    # Calculate small-world properties
    # Compare with random graph of same size and edge count
    n = G.number_of_nodes()
    m = G.number_of_edges()
    
    if n > 1 and m > 0:
        # Create random graph with same density
        random_graph = nx.erdos_renyi_graph(n, 2*m/(n*(n-1)))
        
        if nx.is_connected(random_graph):
            random_avg_path = nx.average_shortest_path_length(random_graph)
            random_clustering = nx.average_clustering(random_graph)
        else:
            if random_graph.number_of_nodes() > 0:
                largest_random_cc = max(nx.connected_components(random_graph), key=len)
                random_subgraph = random_graph.subgraph(largest_random_cc)
                random_avg_path = nx.average_shortest_path_length(random_subgraph) if random_subgraph.number_of_edges() > 0 else 0
            else:
                random_avg_path = 0
            random_clustering = nx.average_clustering(random_graph)
    else:
        random_avg_path = 0
        random_clustering = 0
    
    actual_clustering = nx.average_clustering(G)
    small_world_sigma = (actual_clustering / random_clustering) / (avg_path_length / random_avg_path) if random_clustering > 0 and random_avg_path > 0 else 0

    # Create plots
    organism_suffix = "_excl_synthetic" if exclude_synthetic else "_incl_synthetic"
    organism_title = "Excluding Synthetic Organisms" if exclude_synthetic else "Including Synthetic Organisms"
    
    plt.figure(figsize=(15, 5))
    
    # Subplot 1: Path length distribution
    plt.subplot(1, 3, 1)
    if path_length_counts:
        lengths = sorted(path_length_counts.keys())
        counts = [path_length_counts[l] for l in lengths]
        plt.bar(lengths, counts, alpha=0.7)
    plt.title(f"Protein Network: Shortest Path Length Distribution\n(Single Mutation Steps) - {organism_title}")
    plt.xlabel("Path Length")
    plt.ylabel("Frequency")

    # Subplot 2: Path length vs frequency (log scale)
    plt.subplot(1, 3, 2)
    if path_length_counts:
        lengths = sorted(path_length_counts.keys())
        counts = [path_length_counts[l] for l in lengths]
        if counts and max(counts) > 0:
            plt.semilogy(lengths, counts, 'bo-')
        plt.grid(True, alpha=0.3)
    plt.title(f"Protein Network: Path Length Distribution (Log Scale)\n(Single Mutation Steps) - {organism_title}")
    plt.xlabel("Path Length")
    plt.ylabel("Frequency (log)")

    # Subplot 3: Distance statistics
    plt.subplot(1, 3, 3)
    distance_stats = {
        'Diameter': diameter,
        'Avg Path\nLength': avg_path_length,
        'Efficiency': efficiency,
        'Small World\nSigma': small_world_sigma
    }
    plt.bar(distance_stats.keys(), distance_stats.values(), color='lightblue')
    plt.title(f"Protein Network: Distance Statistics\n(Single Mutation Steps) - {organism_title}")
    plt.ylabel("Value")
    plt.xticks(rotation=45)

    # Bottom row plots removed as requested

    plt.tight_layout()
    plt.savefig(f"{output_path}/long_distance_interference_analysis{organism_suffix}.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Save distance analysis statistics
    with open(f"{output_path}/long_distance_interference_statistics{organism_suffix}.txt", 'w') as f:
        f.write("PROTEIN NETWORK: LONG DISTANCE INTERFERENCE ANALYSIS\n")
        f.write("===================================================\n\n")
        f.write(f"ORGANISM INCLUSION: {organism_title}\n\n")
        f.write("EDGE CRITERIA:\n")
        f.write("Edges represent proteins connected by single mutations only\n")
        f.write("(proteins that differ by exactly one amino acid position)\n\n")
        f.write("PATH LENGTH INTERPRETATION:\n")
        f.write("Path length = minimum number of single mutations required\n")
        f.write("to connect two proteins via evolutionary intermediates.\n\n")
        f.write(f"Network connectivity: {'Connected' if nx.is_connected(G) else 'Disconnected'}\n")
        f.write(f"Number of connected components: {nx.number_connected_components(G)}\n")
        if G.number_of_nodes() > 0:
            f.write(f"Largest component size: {len(max(nx.connected_components(G), key=len))}\n\n")
        else:
            f.write("Largest component size: 0\n\n")
        
        f.write(f"DISTANCE STATISTICS:\n")
        f.write(f"Diameter: {diameter}\n")
        f.write(f"Average shortest path length: {avg_path_length:.3f}\n")
        f.write(f"Network efficiency: {efficiency:.6f}\n\n")
        
        f.write(f"SMALL-WORLD PROPERTIES:\n")
        f.write(f"Actual clustering coefficient: {actual_clustering:.6f}\n")
        f.write(f"Random clustering coefficient: {random_clustering:.6f}\n")
        f.write(f"Actual average path length: {avg_path_length:.3f}\n")
        f.write(f"Random average path length: {random_avg_path:.3f}\n")
        f.write(f"Small-world sigma: {small_world_sigma:.3f}\n")
        
        if small_world_sigma > 1:
            f.write("Network exhibits small-world properties\n")
        else:
            f.write("Network does not exhibit small-world properties\n")

    LOGGER.info(f"Long distance analysis completed ({'excl' if exclude_synthetic else 'incl'} synthetic). Diameter: {diameter}, Avg path length: {avg_path_length:.3f}")
    return {
        'diameter': diameter,
        'avg_path_length': avg_path_length,
        'efficiency': efficiency,
        'small_world_sigma': small_world_sigma
    }
