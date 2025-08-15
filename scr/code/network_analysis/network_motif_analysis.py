import logging
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx


def analyze_graph_motifs(G, output_path, exclude_synthetic=True, suffix_override=None):
    """Search for and analyze graph motifs
    
    Note: Analyzes motifs in single-mutation protein network.
    Triangles = groups of 3 proteins all within single mutations of each other.
    Stars = hub proteins with many single-mutation neighbors.
    
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
    
    # Count different motifs
    motif_counts = {}
    
    # Triangles (3-cliques)
    triangles = list(nx.enumerate_all_cliques(G))
    triangle_count = len([clique for clique in triangles if len(clique) == 3])
    motif_counts['Triangles'] = triangle_count
    
    # Stars (nodes with degree > threshold connected to leaves)
    star_count = 0
    star_sizes = []
    degree_threshold = 5  # Consider nodes with degree > 5 as potential star centers
    
    for node in G.nodes():
        degree = G.degree(node)
        if degree > degree_threshold:
            neighbors = list(G.neighbors(node))
            # Check if this forms a star (neighbors should have low connectivity among themselves)
            neighbor_edges = 0
            for i in range(len(neighbors)):
                for j in range(i+1, len(neighbors)):
                    if G.has_edge(neighbors[i], neighbors[j]):
                        neighbor_edges += 1
            
            # If neighbors have few connections among themselves, it's star-like
            max_possible_edges = len(neighbors) * (len(neighbors) - 1) / 2
            if max_possible_edges > 0 and neighbor_edges / max_possible_edges < 0.3:
                star_count += 1
                star_sizes.append(degree)
    
    motif_counts['Stars'] = star_count
    
    # Paths of length 3 (simple paths)
    path_3_count = 0
    for node in G.nodes():
        neighbors = list(G.neighbors(node))
        for neighbor in neighbors:
            second_neighbors = list(G.neighbors(neighbor))
            for second_neighbor in second_neighbors:
                if second_neighbor != node and not G.has_edge(node, second_neighbor):
                    path_3_count += 1
    
    path_3_count = path_3_count // 2  # Each path is counted twice
    motif_counts['Paths (length 3)'] = path_3_count
    
    # Squares (4-cycles)
    squares = [cycle for cycle in nx.cycle_basis(G) if len(cycle) == 4]
    motif_counts['Squares'] = len(squares)
    
    # Create motif analysis plots
    if suffix_override:
        organism_suffix = f"_{suffix_override}"
        organism_title = "Heterogeneous Network" if "heterogeneous" in suffix_override else ("Excluding Synthetic Organisms" if exclude_synthetic else "Including Synthetic Organisms")
    else:
        organism_suffix = "_excl_synthetic" if exclude_synthetic else "_incl_synthetic"
        organism_title = "Excluding Synthetic Organisms" if exclude_synthetic else "Including Synthetic Organisms"
    
    plt.figure(figsize=(12, 8))
    
    # Subplot 1: Motif counts
    plt.subplot(2, 2, 1)
    if motif_counts:
        plt.bar(motif_counts.keys(), motif_counts.values(), color='skyblue')
        for i, v in enumerate(motif_counts.values()):
            plt.text(i, v + max(motif_counts.values()) * 0.01, str(v), ha='center', va='bottom')
    plt.title(f"Protein Network: Graph Motif Counts\n(Single Mutation Edges) - {organism_title}")
    plt.ylabel("Count")
    plt.xticks(rotation=45)

    # Subplot 2: Star size distribution
    plt.subplot(2, 2, 2)
    if star_sizes:
        plt.hist(star_sizes, bins=20, alpha=0.7, edgecolor='black')
        plt.xlabel("Star Size (Degree)")
        plt.ylabel("Frequency")
    else:
        plt.text(0.5, 0.5, "No stars found", ha='center', va='center', transform=plt.gca().transAxes)
    plt.title(f"Protein Network: Star Size Distribution\n(Single Mutation Edges) - {organism_title}")

    # Subplot 3: Clustering coefficient distribution
    plt.subplot(2, 2, 3)
    clustering_coeffs = list(nx.clustering(G).values())
    if clustering_coeffs:
        plt.hist(clustering_coeffs, bins=30, alpha=0.7, edgecolor='black')
    plt.title(f"Protein Network: Clustering Coefficient Distribution\n(Single Mutation Edges) - {organism_title}")
    plt.xlabel("Clustering Coefficient")
    plt.ylabel("Frequency")

    # Subplot 4: Motif density comparison
    plt.subplot(2, 2, 4)
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    
    # Calculate theoretical maximums for normalization
    max_triangles = (n_nodes * (n_nodes - 1) * (n_nodes - 2)) / 6 if n_nodes >= 3 else 0
    max_paths = n_nodes * (n_nodes - 1) * (n_nodes - 2) if n_nodes >= 3 else 0
    
    normalized_motifs = {
        'Triangle density': triangle_count / max_triangles if max_triangles > 0 else 0,
        'Edge density': (2 * n_edges) / (n_nodes * (n_nodes - 1)) if n_nodes > 1 else 0,
        'Clustering coeff': np.mean(clustering_coeffs) if clustering_coeffs else 0
    }
    
    if normalized_motifs:
        plt.bar(normalized_motifs.keys(), normalized_motifs.values(), color='lightgreen')
        plt.xticks(rotation=45)
    plt.title(f"Protein Network: Normalized Motif Densities\n(Single Mutation Edges) - {organism_title}")
    plt.ylabel("Density")

    plt.tight_layout()
    plt.savefig(f"{output_path}/motif_analysis{organism_suffix}.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Save motif statistics
    with open(f"{output_path}/motif_statistics{organism_suffix}.txt", 'w') as f:
        f.write("PROTEIN NETWORK: GRAPH MOTIF ANALYSIS\n")
        f.write("====================================\n\n")
        f.write(f"ORGANISM INCLUSION: {organism_title}\n\n")
        f.write("EDGE CRITERIA:\n")
        f.write("Edges represent proteins connected by single mutations only\n")
        f.write("(proteins that differ by exactly one amino acid position)\n\n")
        f.write("MOTIF INTERPRETATIONS:\n")
        f.write("- Triangles: Groups of 3 proteins all within single mutations\n")
        f.write("- Stars: Hub proteins with many single-mutation neighbors\n")
        f.write("- Paths: Chains of single-mutation connections\n\n")
        f.write(f"Network size: {n_nodes} nodes, {n_edges} edges\n\n")
        
        f.write("MOTIF COUNTS:\n")
        for motif, count in motif_counts.items():
            f.write(f"{motif}: {count}\n")
        
        f.write(f"\nCLUSTERING STATISTICS:\n")
        if clustering_coeffs:
            f.write(f"Mean clustering coefficient: {np.mean(clustering_coeffs):.6f}\n")
            f.write(f"Global clustering coefficient: {nx.transitivity(G):.6f}\n")
        else:
            f.write("No clustering data available\n")
        
        if star_sizes:
            f.write(f"\nSTAR STATISTICS:\n")
            f.write(f"Mean star size: {np.mean(star_sizes):.2f}\n")
            f.write(f"Max star size: {max(star_sizes)}\n")
            f.write(f"Min star size: {min(star_sizes)}\n")

    LOGGER.info(f"Motif analysis completed ({'excl' if exclude_synthetic else 'incl'} synthetic). Found {triangle_count} triangles, {star_count} stars")
    return motif_counts
