import ast
import logging
import os
from collections import Counter, defaultdict
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.embedding_analysis import EmbeddingTool
from pyeed.analysis.mutation_detection import MutationDetection
from pyeed.analysis.standard_numbering import StandardNumberingTool
from scipy import optimize
import seaborn as sns
from pathlib import Path


def setup_logging():
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )
    return logging.getLogger(__name__)


def setup_database_connection():
    """Initialize database connection using environment variables"""
    load_dotenv()
    password = os.getenv("NEO4J_NIKLAS_TEM_THREE")
    if password is None:
        raise ValueError("NEO4J_NIKLAS_TEM_THREE is not set in the .env file.")
    
    uri = "bolt://129.69.129.130:2137"
    user = "neo4j"
    eedb = Pyeed(uri, user=user, password=password)
    return eedb


def create_output_directory(output_path="/home/nab/Niklas/TEM-lactamase/data/001_results/011_Staan_Stats"):
    """Create output directory if it doesn't exist"""
    Path(output_path).mkdir(parents=True, exist_ok=True)
    return output_path


def get_tem_family_data(eedb):
    """Get TEM beta-lactamase family protein data"""
    data_tem_ids = {}
    base_url_tem_family_card = 'http://purl.obolibrary.org/obo/ARO_3000014'

    query = f"""
    MATCH (o:OntologyObject {{name: '{base_url_tem_family_card}'}})-[*1..1]-(n) RETURN n
    """
    
    result_ontology_object = eedb.db.execute_read(query)
    
    for single_tem in result_ontology_object:
        if single_tem['n']['name'] == 'http://purl.obolibrary.org/obo/ARO_3000078':
            continue
        tem_name = single_tem['n']['label']
        tem_url = single_tem['n']['name']

        query_tem_url = f"""
        MATCH (o:OntologyObject {{name: '{tem_url}'}})-[*1..1]-(n:Protein) RETURN n
        """

        result_tem_url = eedb.db.execute_read(query_tem_url)
        if len(result_tem_url) == 0:
            continue
        result_tem_url = result_tem_url[0]

        if 'IdenticalIds' in result_tem_url['n']:
            data_tem_ids[tem_name] = result_tem_url['n']['IdenticalIds'] + [result_tem_url['n']['accession_id']]
        else:
            data_tem_ids[tem_name] = [result_tem_url['n']['accession_id']]
    
    return data_tem_ids


def get_protein_network_data(eedb, name_of_standard_numbering_tool="standard_numbering_pairwise_flagged_proteins"):
    """Query Neo4j to get protein mutation network data"""
    LOGGER = logging.getLogger(__name__)
    
    query = """
    MATCH (p1:Protein)-[r:MUTATION]->(p2:Protein)
    WHERE p1.seq_length = 286 
    AND p2.seq_length = 286
    AND size(r.from_positions) = 1
    WITH p1, p2, r
    MATCH (o1:Organism)-[:ORIGINATES_FROM]-(p1)
    MATCH (o2:Organism)-[:ORIGINATES_FROM]-(p2)
    WHERE o1.taxonomy_id <> 32630
    AND o2.taxonomy_id <> 32630
    WITH p1, p2, r, o1, o2
    MATCH (s:StandardNumbering {name: $name_of_standard_numbering_tool})-[s1:HAS_STANDARD_NUMBERING]-(p1)
    MATCH (s)-[s2:HAS_STANDARD_NUMBERING]-(p2)
    WITH p1, p2, r, s1, s2, o1, o2, rand() AS random_value, elementId(r) AS id_r
    ORDER BY random_value
    RETURN p1.accession_id AS source, 
           p2.accession_id AS target, 
           r.from_positions[0] AS from_pos, 
           r.to_positions[0] AS to_pos, 
           id_r AS id_r, 
           r.from_monomers AS from_monomer, 
           r.to_monomers AS to_monomer, 
           s1.positions AS position_s1, 
           s2.positions AS position_s2,
           o1.name AS source_organism,
           o2.name AS target_organism
    """

    results = eedb.db.execute_read(query, parameters={"name_of_standard_numbering_tool": name_of_standard_numbering_tool})
    LOGGER.info(f"Number of mutation relationships: {len(results)}")
    return results


def create_networkx_graph(results):
    """Create NetworkX graph from mutation data
    
    Creates a protein network where:
    - Nodes: Individual proteins (identified by accession IDs)
    - Edges: Single mutations only (proteins differing by exactly one amino acid)
    
    This represents the single-step mutational landscape of the protein family.
    """
    G = nx.Graph()
    
    for record in results:
        source = record["source"]
        target = record["target"]
        from_pos = record["from_pos"]
        to_pos = record["to_pos"]

        if source not in G:
            G.add_node(source)
        if target not in G:
            G.add_node(target)

        G.add_edge(source, target, from_pos=from_pos, to_pos=to_pos)
    
    return G


def analyze_node_degrees(G, output_path):
    """Perform node degree analysis and create visualizations
    
    Note: Edges represent single mutations only (proteins that differ
    by exactly one amino acid position). Degree = number of proteins
    reachable by single mutations.
    """
    LOGGER = logging.getLogger(__name__)
    
    # Calculate degree sequence and distribution
    degree_sequence = [d for _, d in G.degree()]
    degree_counts = Counter(degree_sequence)
    degrees = sorted(degree_counts.keys())
    counts = [degree_counts[d] for d in degrees]

    # Log-binning for better visualization
    max_degree = max(degrees)
    bin_edges = np.unique(2 ** np.arange(0, int(np.log2(max_degree))+1))
    
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



    # Create plots
    plt.figure(figsize=(12, 8))
    
    # Subplot 1: Regular degree distribution
    plt.subplot(2, 2, 1)
    plt.bar(degrees, counts, alpha=0.7)
    plt.title("Protein Network: Degree Distribution\n(Single Mutation Edges)")
    plt.xlabel("Degree")
    plt.ylabel("Frequency")
    plt.grid(True, alpha=0.3)

    # Subplot 2: Log-log degree distribution with binning
    plt.subplot(2, 2, 2)
    plt.loglog(degrees, counts, "bo-", linewidth=1, markersize=4, alpha=0.6, label='Raw data')
    if len(binned_degrees) > 0:
        plt.loglog(binned_degrees, binned_counts, "rs-", linewidth=2, markersize=8, label='Log-binned')
    plt.grid(True, which="both", ls="--", alpha=0.4)
    plt.title("Protein Network: Degree Distribution (Log-Log Scale)\n(Single Mutation Edges)")
    plt.xlabel("Degree (log)")
    plt.ylabel("Frequency (log)")
    plt.legend()

    # Subplot 3: Cumulative degree distribution
    plt.subplot(2, 2, 3)
    sorted_degrees = sorted(degree_sequence, reverse=True)
    cumulative = np.arange(1, len(sorted_degrees) + 1) / len(sorted_degrees)
    plt.loglog(sorted_degrees, cumulative, "go-", linewidth=2, markersize=4)
    plt.title("Protein Network: Cumulative Degree Distribution\n(Single Mutation Edges)")
    plt.xlabel("Degree (log)")
    plt.ylabel("P(X ≥ k)")
    plt.grid(True, which="both", ls="--", alpha=0.4)

    # Subplot 4: Degree statistics
    plt.subplot(2, 2, 4)
    degree_stats = {
        'Mean': np.mean(degree_sequence),
        'Median': np.median(degree_sequence),
        'Std': np.std(degree_sequence),
        'Max': max(degree_sequence),
        'Min': min(degree_sequence)
    }
    plt.bar(degree_stats.keys(), degree_stats.values(), color='lightblue')
    plt.title("Protein Network: Degree Statistics\n(Single Mutation Edges)")
    plt.ylabel("Value")
    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig(f"{output_path}/degree_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Save degree statistics to file
    with open(f"{output_path}/degree_statistics.txt", 'w') as f:
        f.write("PROTEIN NETWORK: NODE DEGREE ANALYSIS\n")
        f.write("====================================\n\n")
        f.write("EDGE CRITERIA:\n")
        f.write("Edges represent proteins connected by single mutations only\n")
        f.write("(proteins that differ by exactly one amino acid position)\n\n")
        f.write("INTERPRETATION:\n")
        f.write("Node degree = number of proteins reachable by single mutations\n")
        f.write("High degree = protein with many single-step mutational neighbors\n\n")
        f.write(f"Number of nodes: {G.number_of_nodes()}\n")
        f.write(f"Number of edges: {G.number_of_edges()}\n")
        f.write(f"Mean degree: {np.mean(degree_sequence):.3f}\n")
        f.write(f"Median degree: {np.median(degree_sequence):.3f}\n")
        f.write(f"Standard deviation: {np.std(degree_sequence):.3f}\n")
        f.write(f"Maximum degree: {max(degree_sequence)}\n")
        f.write(f"Minimum degree: {min(degree_sequence)}\n")
        f.write(f"\nLOG BINNING INFORMATION:\n")
        f.write(f"Number of log bins: {len(binned_degrees)}\n")
        f.write(f"Bin edges: {bin_edges.tolist() if len(bin_edges) > 0 else 'None'}\n")
        f.write(f"Binned degrees: {binned_degrees.tolist() if len(binned_degrees) > 0 else 'None'}\n")
        f.write(f"Binned counts: {binned_counts.tolist() if len(binned_counts) > 0 else 'None'}\n")

    LOGGER.info(f"Degree analysis completed. Mean degree: {np.mean(degree_sequence):.3f}")
    return degree_stats


def analyze_betweenness_centrality(G, output_path, data_tem_ids):
    """Calculate and analyze betweenness centrality
    
    Note: Edges in this protein network represent single mutations only
    (proteins that differ by exactly one amino acid position).
    Betweenness centrality measures how important each protein is as an
    intermediary in single-step mutation pathways.
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
    plt.figure(figsize=(15, 10))
    
    # Subplot 1: Top betweenness centrality nodes
    plt.subplot(2, 2, 1)
    plt.bar(range(len(labeled_nodes)), labeled_nodes.values())
    plt.title("Protein Network: Top 20 Nodes by Betweenness Centrality\n(Single Mutation Edges)")
    plt.xlabel("Node")
    plt.ylabel("Betweenness Centrality")
    plt.xticks(range(len(labeled_nodes)), labeled_nodes.keys(), rotation=45, ha='right')

    # Subplot 2: Betweenness centrality distribution
    plt.subplot(2, 2, 2)
    centrality_values = list(betweenness_centrality.values())
    plt.hist(centrality_values, bins=50, alpha=0.7, edgecolor='black')
    plt.title("Protein Network: Betweenness Centrality Distribution\n(Single Mutation Edges)")
    plt.xlabel("Betweenness Centrality")
    plt.ylabel("Frequency")
    plt.yscale('log')

    # Subplot 3: Centrality vs Degree correlation
    plt.subplot(2, 2, 3)
    degrees = [G.degree(node) for node in G.nodes()]
    centralities = [betweenness_centrality[node] for node in G.nodes()]
    plt.scatter(degrees, centralities, alpha=0.6)
    plt.xlabel("Node Degree")
    plt.ylabel("Betweenness Centrality")
    plt.title("Protein Network: Degree vs Betweenness Centrality\n(Single Mutation Edges)")
    
    # Calculate correlation
    correlation = np.corrcoef(degrees, centralities)[0, 1]
    plt.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
             transform=plt.gca().transAxes, fontsize=12, 
             bbox=dict(boxstyle="round", facecolor='wheat', alpha=0.8))

    # Subplot 4: Centrality statistics
    plt.subplot(2, 2, 4)
    centrality_stats = {
        'Mean': np.mean(centrality_values),
        'Median': np.median(centrality_values),
        'Std': np.std(centrality_values),
        'Max': max(centrality_values),
        '90th percentile': np.percentile(centrality_values, 90)
    }
    plt.bar(centrality_stats.keys(), centrality_stats.values(), color='lightcoral')
    plt.title("Protein Network: Betweenness Centrality Statistics\n(Single Mutation Edges)")
    plt.ylabel("Value")
    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig(f"{output_path}/betweenness_centrality_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Save centrality statistics
    with open(f"{output_path}/betweenness_centrality_statistics.txt", 'w') as f:
        f.write("PROTEIN NETWORK: BETWEENNESS CENTRALITY ANALYSIS\n")
        f.write("===============================================\n\n")
        f.write("EDGE CRITERIA:\n")
        f.write("Edges represent proteins connected by single mutations only\n")
        f.write("(proteins that differ by exactly one amino acid position)\n\n")
        f.write("INTERPRETATION:\n")
        f.write("Betweenness centrality measures how important each protein is\n")
        f.write("as an intermediary in single-step mutation pathways.\n")
        f.write("High centrality = key stepping stone in evolutionary paths.\n\n")
        f.write(f"Mean betweenness centrality: {np.mean(centrality_values):.6f}\n")
        f.write(f"Median betweenness centrality: {np.median(centrality_values):.6f}\n")
        f.write(f"Standard deviation: {np.std(centrality_values):.6f}\n")
        f.write(f"Maximum centrality: {max(centrality_values):.6f}\n")
        f.write(f"Correlation with degree: {correlation:.3f}\n\n")
        
        f.write("TOP 20 NODES BY BETWEENNESS CENTRALITY:\n")
        for i, (node, centrality) in enumerate(top_nodes):
            name = id_to_name.get(node, "Unknown")
            f.write(f"{i+1:2d}. {node} → {name}: {centrality:.6f}\n")

    LOGGER.info(f"Betweenness centrality analysis completed. Max centrality: {max(centrality_values):.6f}")
    return betweenness_centrality


def analyze_graph_motifs(G, output_path):
    """Search for and analyze graph motifs
    
    Note: Analyzes motifs in single-mutation protein network.
    Triangles = groups of 3 proteins all within single mutations of each other.
    Stars = hub proteins with many single-mutation neighbors.
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
    plt.figure(figsize=(12, 8))
    
    # Subplot 1: Motif counts
    plt.subplot(2, 2, 1)
    plt.bar(motif_counts.keys(), motif_counts.values(), color='skyblue')
    plt.title("Protein Network: Graph Motif Counts\n(Single Mutation Edges)")
    plt.ylabel("Count")
    plt.xticks(rotation=45)
    for i, v in enumerate(motif_counts.values()):
        plt.text(i, v + max(motif_counts.values()) * 0.01, str(v), ha='center', va='bottom')

    # Subplot 2: Star size distribution
    plt.subplot(2, 2, 2)
    if star_sizes:
        plt.hist(star_sizes, bins=20, alpha=0.7, edgecolor='black')
        plt.title("Protein Network: Star Size Distribution\n(Single Mutation Edges)")
        plt.xlabel("Star Size (Degree)")
        plt.ylabel("Frequency")
    else:
        plt.text(0.5, 0.5, "No stars found", ha='center', va='center', transform=plt.gca().transAxes)
        plt.title("Protein Network: Star Size Distribution\n(Single Mutation Edges)")

    # Subplot 3: Clustering coefficient distribution
    plt.subplot(2, 2, 3)
    clustering_coeffs = list(nx.clustering(G).values())
    plt.hist(clustering_coeffs, bins=30, alpha=0.7, edgecolor='black')
    plt.title("Protein Network: Clustering Coefficient Distribution\n(Single Mutation Edges)")
    plt.xlabel("Clustering Coefficient")
    plt.ylabel("Frequency")

    # Subplot 4: Motif density comparison
    plt.subplot(2, 2, 4)
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    
    # Calculate theoretical maximums for normalization
    max_triangles = (n_nodes * (n_nodes - 1) * (n_nodes - 2)) / 6
    max_paths = n_nodes * (n_nodes - 1) * (n_nodes - 2)
    
    normalized_motifs = {
        'Triangle density': triangle_count / max_triangles if max_triangles > 0 else 0,
        'Edge density': (2 * n_edges) / (n_nodes * (n_nodes - 1)) if n_nodes > 1 else 0,
        'Clustering coeff': np.mean(clustering_coeffs) if clustering_coeffs else 0
    }
    
    plt.bar(normalized_motifs.keys(), normalized_motifs.values(), color='lightgreen')
    plt.title("Protein Network: Normalized Motif Densities\n(Single Mutation Edges)")
    plt.ylabel("Density")
    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig(f"{output_path}/motif_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Save motif statistics
    with open(f"{output_path}/motif_statistics.txt", 'w') as f:
        f.write("PROTEIN NETWORK: GRAPH MOTIF ANALYSIS\n")
        f.write("====================================\n\n")
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
        f.write(f"Mean clustering coefficient: {np.mean(clustering_coeffs):.6f}\n")
        f.write(f"Global clustering coefficient: {nx.transitivity(G):.6f}\n")
        
        if star_sizes:
            f.write(f"\nSTAR STATISTICS:\n")
            f.write(f"Mean star size: {np.mean(star_sizes):.2f}\n")
            f.write(f"Max star size: {max(star_sizes)}\n")
            f.write(f"Min star size: {min(star_sizes)}\n")

    LOGGER.info(f"Motif analysis completed. Found {triangle_count} triangles, {star_count} stars")
    return motif_counts


def analyze_long_distance_interference(G, output_path):
    """Investigate long distance interference in the graph
    
    Note: Analyzes path lengths in single-mutation protein network.
    Path length = minimum number of single mutations required to connect two proteins.
    This reveals the mutational distance and evolutionary relationships.
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
        
        diameter = max(path_lengths)
        avg_path_length = np.mean(path_lengths)
        
    else:
        # For disconnected graph, analyze largest component
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
    
    # Calculate efficiency and resistance distance approximation
    efficiency = 1 / avg_path_length if avg_path_length > 0 else 0
    
    # Analyze path length distribution
    path_length_counts = Counter(path_lengths)
    
    # Calculate small-world properties
    # Compare with random graph of same size and edge count
    n = G.number_of_nodes()
    m = G.number_of_edges()
    random_graph = nx.erdos_renyi_graph(n, 2*m/(n*(n-1)))
    
    if nx.is_connected(random_graph):
        random_avg_path = nx.average_shortest_path_length(random_graph)
        random_clustering = nx.average_clustering(random_graph)
    else:
        largest_random_cc = max(nx.connected_components(random_graph), key=len)
        random_subgraph = random_graph.subgraph(largest_random_cc)
        random_avg_path = nx.average_shortest_path_length(random_subgraph)
        random_clustering = nx.average_clustering(random_graph)
    
    actual_clustering = nx.average_clustering(G)
    small_world_sigma = (actual_clustering / random_clustering) / (avg_path_length / random_avg_path) if random_clustering > 0 and random_avg_path > 0 else 0

    # Create plots
    plt.figure(figsize=(15, 10))
    
    # Subplot 1: Path length distribution
    plt.subplot(2, 3, 1)
    lengths = sorted(path_length_counts.keys())
    counts = [path_length_counts[l] for l in lengths]
    plt.bar(lengths, counts, alpha=0.7)
    plt.title("Protein Network: Shortest Path Length Distribution\n(Single Mutation Steps)")
    plt.xlabel("Path Length")
    plt.ylabel("Frequency")

    # Subplot 2: Path length vs frequency (log scale)
    plt.subplot(2, 3, 2)
    plt.semilogy(lengths, counts, 'bo-')
    plt.title("Protein Network: Path Length Distribution (Log Scale)\n(Single Mutation Steps)")
    plt.xlabel("Path Length")
    plt.ylabel("Frequency (log)")
    plt.grid(True, alpha=0.3)

    # Subplot 3: Distance statistics
    plt.subplot(2, 3, 3)
    distance_stats = {
        'Diameter': diameter,
        'Avg Path\nLength': avg_path_length,
        'Efficiency': efficiency,
        'Small World\nSigma': small_world_sigma
    }
    plt.bar(distance_stats.keys(), distance_stats.values(), color='lightblue')
    plt.title("Protein Network: Distance Statistics\n(Single Mutation Steps)")
    plt.ylabel("Value")
    plt.xticks(rotation=45)

    # Subplot 4: Component analysis
    plt.subplot(2, 3, 4)
    components = list(nx.connected_components(G))
    component_sizes = [len(comp) for comp in components]
    plt.hist(component_sizes, bins=20, alpha=0.7, edgecolor='black')
    plt.title("Protein Network: Connected Component Sizes\n(Single Mutation Edges)")
    plt.xlabel("Component Size")
    plt.ylabel("Frequency")
    plt.yscale('log')

    # Subplot 5: Resistance distance visualization (sample)
    plt.subplot(2, 3, 5)
    # Sample a subset of nodes for resistance distance calculation
    sample_nodes = list(G.nodes())[:min(50, len(G.nodes()))]
    sample_subgraph = G.subgraph(sample_nodes)
    
    if nx.is_connected(sample_subgraph):
        # Calculate resistance distances for small sample
        resistance_distances = []
        for i, node1 in enumerate(sample_nodes):
            for node2 in sample_nodes[i+1:]:
                try:
                    path_length = nx.shortest_path_length(sample_subgraph, node1, node2)
                    resistance_distances.append(path_length)
                except nx.NetworkXNoPath:
                    pass
        
        if resistance_distances:
            plt.hist(resistance_distances, bins=15, alpha=0.7, edgecolor='black')
            plt.title(f"Protein Network: Sample Resistance Distances\n(Single Mutation Steps, n={len(sample_nodes)})")
            plt.xlabel("Resistance Distance")
            plt.ylabel("Frequency")
    else:
        plt.text(0.5, 0.5, "Sample not connected", ha='center', va='center', 
                transform=plt.gca().transAxes)
        plt.title("Protein Network: Sample Resistance Distances\n(Single Mutation Steps)")

    # Subplot 6: Small-world comparison
    plt.subplot(2, 3, 6)
    comparison_data = {
        'Actual\nClustering': actual_clustering,
        'Random\nClustering': random_clustering,
        'Actual\nPath Length': avg_path_length,
        'Random\nPath Length': random_avg_path
    }
    
    colors = ['blue', 'lightblue', 'red', 'lightcoral']
    plt.bar(comparison_data.keys(), comparison_data.values(), color=colors)
    plt.title("Protein Network: Small-World Comparison\n(Single Mutation Network)")
    plt.ylabel("Value")
    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig(f"{output_path}/long_distance_interference_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Save distance analysis statistics
    with open(f"{output_path}/long_distance_interference_statistics.txt", 'w') as f:
        f.write("PROTEIN NETWORK: LONG DISTANCE INTERFERENCE ANALYSIS\n")
        f.write("===================================================\n\n")
        f.write("EDGE CRITERIA:\n")
        f.write("Edges represent proteins connected by single mutations only\n")
        f.write("(proteins that differ by exactly one amino acid position)\n\n")
        f.write("PATH LENGTH INTERPRETATION:\n")
        f.write("Path length = minimum number of single mutations required\n")
        f.write("to connect two proteins via evolutionary intermediates.\n\n")
        f.write(f"Network connectivity: {'Connected' if nx.is_connected(G) else 'Disconnected'}\n")
        f.write(f"Number of connected components: {nx.number_connected_components(G)}\n")
        f.write(f"Largest component size: {len(max(nx.connected_components(G), key=len))}\n\n")
        
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

    LOGGER.info(f"Long distance analysis completed. Diameter: {diameter}, Avg path length: {avg_path_length:.3f}")
    return {
        'diameter': diameter,
        'avg_path_length': avg_path_length,
        'efficiency': efficiency,
        'small_world_sigma': small_world_sigma
    }


def analyze_mutation_lengths_by_phenotype(eedb, results, output_path):
    """Analyze mutation lengths for each phenotype class"""
    LOGGER = logging.getLogger(__name__)
    
    # Load phenotype data
    try:
        df = pd.read_csv("/home/nab/Niklas/TEM-lactamase/data/002_combined_data/TEM_lactamase_with_dna_accession_id.csv")
        phenotype_dict = dict(zip(df['protein_id_database'].dropna(), df['phenotype'].dropna()))
    except FileNotFoundError:
        LOGGER.warning("Phenotype data file not found. Creating mock data.")
        # Create mock phenotype data based on available results
        unique_proteins = set([r['source'] for r in results] + [r['target'] for r in results])
        phenotypes = ['2b', '2be', '2br', '2ber', 'unknown']
        phenotype_dict = {protein: np.random.choice(phenotypes) for protein in unique_proteins}

    # Analyze mutation lengths by phenotype
    phenotype_mutations = defaultdict(list)
    mutation_position_analysis = defaultdict(list)
    
    for result in results:
        source = result['source']
        target = result['target']
        from_pos = result['from_pos']
        to_pos = result['to_pos']
        
        # Get phenotypes
        source_phenotype = phenotype_dict.get(source, 'unknown')
        target_phenotype = phenotype_dict.get(target, 'unknown')
        
        # Store mutation information
        mutation_length = abs(to_pos - from_pos) if from_pos != to_pos else 0
        
        phenotype_mutations[source_phenotype].append(mutation_length)
        phenotype_mutations[target_phenotype].append(mutation_length)
        
        mutation_position_analysis[source_phenotype].append(from_pos)
        mutation_position_analysis[target_phenotype].append(to_pos)

    # Calculate statistics for each phenotype
    phenotype_stats = {}
    for phenotype, lengths in phenotype_mutations.items():
        if lengths:
            phenotype_stats[phenotype] = {
                'count': len(lengths),
                'mean_length': np.mean(lengths),
                'median_length': np.median(lengths),
                'std_length': np.std(lengths),
                'max_length': max(lengths),
                'min_length': min(lengths)
            }

    # Create comprehensive plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1: Mutation length distribution by phenotype
    phenotypes = list(phenotype_stats.keys())
    all_lengths = [phenotype_mutations[p] for p in phenotypes]
    
    axes[0, 0].boxplot(all_lengths, labels=phenotypes)
    axes[0, 0].set_title("Mutation Length Distribution by Phenotype")
    axes[0, 0].set_ylabel("Mutation Length")
    axes[0, 0].tick_params(axis='x', rotation=45)

    # Plot 2: Mean mutation lengths
    mean_lengths = [phenotype_stats[p]['mean_length'] for p in phenotypes]
    colors = plt.cm.Set3(np.linspace(0, 1, len(phenotypes)))
    
    bars = axes[0, 1].bar(phenotypes, mean_lengths, color=colors)
    axes[0, 1].set_title("Mean Mutation Length by Phenotype")
    axes[0, 1].set_ylabel("Mean Mutation Length")
    axes[0, 1].tick_params(axis='x', rotation=45)
    
    # Add value labels on bars
    for bar, value in zip(bars, mean_lengths):
        axes[0, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                       f'{value:.2f}', ha='center', va='bottom')

    # Plot 3: Mutation count by phenotype
    counts = [phenotype_stats[p]['count'] for p in phenotypes]
    axes[0, 2].bar(phenotypes, counts, color=colors)
    axes[0, 2].set_title("Number of Mutations by Phenotype")
    axes[0, 2].set_ylabel("Count")
    axes[0, 2].tick_params(axis='x', rotation=45)

    # Plot 4: Mutation position heatmap
    position_matrix = []
    position_labels = []
    for phenotype in phenotypes:
        positions = mutation_position_analysis[phenotype]
        if positions:
            hist, bin_edges = np.histogram(positions, bins=20, range=(1, 286))
            position_matrix.append(hist)
            position_labels.append(phenotype)
    
    if position_matrix:
        im = axes[1, 0].imshow(position_matrix, cmap='YlOrRd', aspect='auto')
        axes[1, 0].set_title("Mutation Position Heatmap by Phenotype")
        axes[1, 0].set_xlabel("Position Bins")
        axes[1, 0].set_ylabel("Phenotype")
        axes[1, 0].set_yticks(range(len(position_labels)))
        axes[1, 0].set_yticklabels(position_labels)
        plt.colorbar(im, ax=axes[1, 0])

    # Plot 5: Standard deviation comparison
    std_lengths = [phenotype_stats[p]['std_length'] for p in phenotypes]
    axes[1, 1].bar(phenotypes, std_lengths, color=colors, alpha=0.7)
    axes[1, 1].set_title("Mutation Length Variability by Phenotype")
    axes[1, 1].set_ylabel("Standard Deviation")
    axes[1, 1].tick_params(axis='x', rotation=45)

    # Plot 6: Range comparison (max - min)
    ranges = [phenotype_stats[p]['max_length'] - phenotype_stats[p]['min_length'] for p in phenotypes]
    axes[1, 2].bar(phenotypes, ranges, color=colors, alpha=0.7)
    axes[1, 2].set_title("Mutation Length Range by Phenotype")
    axes[1, 2].set_ylabel("Range (Max - Min)")
    axes[1, 2].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(f"{output_path}/mutation_length_phenotype_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Save detailed statistics
    with open(f"{output_path}/mutation_length_phenotype_statistics.txt", 'w') as f:
        f.write("MUTATION LENGTH ANALYSIS BY PHENOTYPE\n")
        f.write("====================================\n\n")
        
        for phenotype, stats in phenotype_stats.items():
            f.write(f"PHENOTYPE: {phenotype}\n")
            f.write(f"  Number of mutations: {stats['count']}\n")
            f.write(f"  Mean length: {stats['mean_length']:.3f}\n")
            f.write(f"  Median length: {stats['median_length']:.3f}\n")
            f.write(f"  Standard deviation: {stats['std_length']:.3f}\n")
            f.write(f"  Range: {stats['min_length']} - {stats['max_length']}\n\n")
        
        # Statistical comparisons
        f.write("STATISTICAL COMPARISONS:\n")
        phenotype_pairs = [(p1, p2) for i, p1 in enumerate(phenotypes) for p2 in phenotypes[i+1:]]
        
        for p1, p2 in phenotype_pairs:
            lengths1 = phenotype_mutations[p1]
            lengths2 = phenotype_mutations[p2]
            if len(lengths1) > 5 and len(lengths2) > 5:
                from scipy.stats import mannwhitneyu
                try:
                    statistic, p_value = mannwhitneyu(lengths1, lengths2, alternative='two-sided')
                    f.write(f"{p1} vs {p2}: Mann-Whitney U test p-value = {p_value:.6f}\n")
                except:
                    f.write(f"{p1} vs {p2}: Statistical test failed\n")

    LOGGER.info(f"Mutation length analysis completed for {len(phenotype_stats)} phenotypes")
    return phenotype_stats


def analyze_database_schema(eedb):
    """Analyze the complete database schema to understand available node types and relationships"""
    LOGGER = logging.getLogger(__name__)
    
    # Get all node labels (types)
    node_labels_query = """
    CALL db.labels() YIELD label
    RETURN label
    """
    node_labels = eedb.db.execute_read(node_labels_query)
    
    # Get all relationship types  
    relationship_types_query = """
    CALL db.relationshipTypes() YIELD relationshipType
    RETURN relationshipType
    """
    relationship_types = eedb.db.execute_read(relationship_types_query)
    
    # Count nodes by type
    node_counts = {}
    for label_record in node_labels:
        label = label_record['label']
        count_query = f"MATCH (n:{label}) RETURN count(n) as count"
        try:
            count_result = eedb.db.execute_read(count_query)
            node_counts[label] = count_result[0]['count'] if count_result else 0
        except:
            node_counts[label] = 0
    
    LOGGER.info(f"Found {len(node_labels)} node types and {len(relationship_types)} relationship types")
    
    return {
        'node_labels': [record['label'] for record in node_labels],
        'relationship_types': [record['relationshipType'] for record in relationship_types],
        'node_counts': node_counts
    }


def get_heterogeneous_network_data(eedb):
    """Query Neo4j to get heterogeneous network data including multiple node types"""
    LOGGER = logging.getLogger(__name__)
    
    # Query for Protein-Protein MUTATION relationships
    protein_mutation_query = """
    MATCH (p1:Protein)-[r:MUTATION]->(p2:Protein)
    WHERE p1.seq_length = 286 AND p2.seq_length = 286
    RETURN 
        'Protein' as source_type, p1.accession_id AS source_id,
        'Protein' as target_type, p2.accession_id AS target_id,
        'MUTATION' as relationship_type,
        r.from_positions as mutation_positions
    """
    
    # Query for Protein-Organism relationships
    protein_organism_query = """
    MATCH (p:Protein)-[r:ORIGINATES_FROM]->(o:Organism)
    WHERE p.seq_length = 286
    RETURN 
        'Protein' as source_type, p.accession_id AS source_id,
        'Organism' as target_type, o.taxonomy_id AS target_id,
        'ORIGINATES_FROM' as relationship_type,
        o.name as organism_name
    """
    
    # Query for DNA-Protein relationships
    dna_protein_query = """
    MATCH (d:DNA)-[r:ENCODES]->(p:Protein)
    WHERE p.seq_length = 286
    RETURN 
        'DNA' as source_type, d.accession_id AS source_id,
        'Protein' as target_type, p.accession_id AS target_id,
        'ENCODES' as relationship_type,
        d.gc_content as gc_content
    """
    
    # Query for Protein-OntologyObject relationships (resistance mechanisms)
    protein_ontology_query = """
    MATCH (p:Protein)-[r]-(o:OntologyObject)
    WHERE p.seq_length = 286
    RETURN 
        'Protein' as source_type, p.accession_id AS source_id,
        'OntologyObject' as target_type, o.name AS target_id,
        type(r) as relationship_type,
        o.label as ontology_label
    LIMIT 1000
    """
    
    # Execute all queries
    all_relationships = []
    
    queries = [
        ("protein_mutations", protein_mutation_query),
        ("protein_organisms", protein_organism_query), 
        ("dna_proteins", dna_protein_query),
        ("protein_ontology", protein_ontology_query)
    ]
    
    for query_name, query in queries:
        try:
            results = eedb.db.execute_read(query)
            LOGGER.info(f"{query_name}: {len(results)} relationships")
            all_relationships.extend(results)
        except Exception as e:
            LOGGER.warning(f"Query {query_name} failed: {e}")
    
    LOGGER.info(f"Total heterogeneous relationships collected: {len(all_relationships)}")
    return all_relationships


def create_heterogeneous_networkx_graph(relationships):
    """Create a heterogeneous NetworkX graph from multi-type relationships"""
    G = nx.Graph()
    
    # Track node types and relationship types
    node_types = {}
    edge_types = {}
    
    for rel in relationships:
        source = f"{rel['source_type']}_{rel['source_id']}"
        target = f"{rel['target_type']}_{rel['target_id']}"
        rel_type = rel['relationship_type']
        
        # Add nodes with type information
        if source not in G:
            G.add_node(source, 
                      node_type=rel['source_type'],
                      original_id=rel['source_id'])
            node_types[source] = rel['source_type']
            
        if target not in G:
            G.add_node(target,
                      node_type=rel['target_type'], 
                      original_id=rel['target_id'])
            node_types[target] = rel['target_type']
        
        # Add edge with relationship type
        G.add_edge(source, target, relationship_type=rel_type)
        edge_types[(source, target)] = rel_type
    
    # Add node type and edge type as graph attributes
    nx.set_node_attributes(G, node_types, 'node_type')
    nx.set_edge_attributes(G, edge_types, 'relationship_type')
    
    return G


def analyze_heterogeneous_network_structure(G_het, output_path):
    """Analyze the structure of the heterogeneous network"""
    LOGGER = logging.getLogger(__name__)
    
    # Get node type distribution
    node_types = nx.get_node_attributes(G_het, 'node_type')
    node_type_counts = Counter(node_types.values())
    
    # Get edge type distribution  
    edge_types = nx.get_edge_attributes(G_het, 'relationship_type')
    edge_type_counts = Counter(edge_types.values())
    
    # Analyze connectivity patterns between node types
    type_connectivity = defaultdict(lambda: defaultdict(int))
    for u, v in G_het.edges():
        source_type = G_het.nodes[u]['node_type']
        target_type = G_het.nodes[v]['node_type']
        type_connectivity[source_type][target_type] += 1
        if source_type != target_type:  # For undirected, count both directions
            type_connectivity[target_type][source_type] += 1
    
    # Create comprehensive visualization
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1: Node type distribution
    axes[0, 0].bar(node_type_counts.keys(), node_type_counts.values(), color='skyblue')
    axes[0, 0].set_title("Heterogeneous: Node Type Distribution")
    axes[0, 0].set_ylabel("Count")
    axes[0, 0].tick_params(axis='x', rotation=45)
    
    # Plot 2: Edge type distribution  
    axes[0, 1].bar(edge_type_counts.keys(), edge_type_counts.values(), color='lightcoral')
    axes[0, 1].set_title("Heterogeneous: Relationship Type Distribution")
    axes[0, 1].set_ylabel("Count")
    axes[0, 1].tick_params(axis='x', rotation=45)
    
    # Plot 3: Network connectivity heatmap
    type_names = list(node_type_counts.keys())
    connectivity_matrix = np.zeros((len(type_names), len(type_names)))
    for i, type1 in enumerate(type_names):
        for j, type2 in enumerate(type_names):
            connectivity_matrix[i, j] = type_connectivity[type1][type2]
    
    im = axes[0, 2].imshow(connectivity_matrix, cmap='YlOrRd')
    axes[0, 2].set_title("Heterogeneous: Inter-Type Connectivity Matrix")
    axes[0, 2].set_xticks(range(len(type_names)))
    axes[0, 2].set_yticks(range(len(type_names)))
    axes[0, 2].set_xticklabels(type_names, rotation=45)
    axes[0, 2].set_yticklabels(type_names)
    plt.colorbar(im, ax=axes[0, 2])
    
    # Plot 4: Degree distribution by node type
    degrees_by_type = defaultdict(list)
    for node in G_het.nodes():
        node_type = G_het.nodes[node]['node_type']
        degree = G_het.degree(node)
        degrees_by_type[node_type].append(degree)
    
    type_colors = plt.cm.Set3(np.linspace(0, 1, len(degrees_by_type)))
    for i, (node_type, degrees) in enumerate(degrees_by_type.items()):
        if degrees:  # Only plot if there are degrees to plot
            axes[1, 0].hist(degrees, alpha=0.7, label=node_type, color=type_colors[i], bins=20)
    axes[1, 0].set_title("Heterogeneous: Degree Distribution by Node Type")
    axes[1, 0].set_xlabel("Degree")
    axes[1, 0].set_ylabel("Frequency")
    axes[1, 0].legend()
    axes[1, 0].set_yscale('log')
    
    # Plot 5: Average degree by node type
    avg_degrees = {node_type: np.mean(degrees) if degrees else 0 for node_type, degrees in degrees_by_type.items()}
    axes[1, 1].bar(avg_degrees.keys(), avg_degrees.values(), color=type_colors)
    axes[1, 1].set_title("Heterogeneous: Average Degree by Node Type")
    axes[1, 1].set_ylabel("Average Degree")
    axes[1, 1].tick_params(axis='x', rotation=45)
    
    # Plot 6: Network density by relationship type
    density_by_rel_type = {}
    for rel_type in edge_type_counts.keys():
        # Create subgraph with only this relationship type
        edges_of_type = [(u, v) for u, v, d in G_het.edges(data=True) 
                        if d.get('relationship_type') == rel_type]
        if edges_of_type:
            subgraph = G_het.edge_subgraph(edges_of_type)
            density_by_rel_type[rel_type] = nx.density(subgraph)
    
    axes[1, 2].bar(density_by_rel_type.keys(), density_by_rel_type.values(), color='lightgreen')
    axes[1, 2].set_title("Heterogeneous: Network Density by Relationship Type")
    axes[1, 2].set_ylabel("Density")
    axes[1, 2].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/heterogeneous_network_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    LOGGER.info("Heterogeneous network structure analysis completed")
    return {
        'node_type_counts': node_type_counts,
        'edge_type_counts': edge_type_counts,
        'avg_degrees': avg_degrees,
        'connectivity_matrix': connectivity_matrix
    }


def compare_network_approaches(G_protein, G_het, output_path, data_tem_ids):
    """Compare protein-only vs heterogeneous network approaches"""
    LOGGER = logging.getLogger(__name__)
    
    # Extract protein subgraph from heterogeneous network
    protein_nodes_het = [n for n in G_het.nodes() if G_het.nodes[n]['node_type'] == 'Protein']
    G_protein_from_het = G_het.subgraph(protein_nodes_het)
    
    # Calculate comparative metrics
    metrics_comparison = {
        'Protein-only Network': {
            'nodes': G_protein.number_of_nodes(),
            'edges': G_protein.number_of_edges(),
            'density': nx.density(G_protein),
            'avg_clustering': nx.average_clustering(G_protein),
            'components': nx.number_connected_components(G_protein),
            'avg_degree': np.mean([d for n, d in G_protein.degree()]) if G_protein.number_of_nodes() > 0 else 0,
            'diameter': 0,  # Will calculate if connected
            'avg_path_length': 0
        },
        'Heterogeneous Network (full)': {
            'nodes': G_het.number_of_nodes(),
            'edges': G_het.number_of_edges(), 
            'density': nx.density(G_het),
            'avg_clustering': nx.average_clustering(G_het),
            'components': nx.number_connected_components(G_het),
            'avg_degree': np.mean([d for n, d in G_het.degree()]) if G_het.number_of_nodes() > 0 else 0,
            'diameter': 0,
            'avg_path_length': 0
        },
        'Protein subgraph from Het.': {
            'nodes': G_protein_from_het.number_of_nodes(),
            'edges': G_protein_from_het.number_of_edges(),
            'density': nx.density(G_protein_from_het),
            'avg_clustering': nx.average_clustering(G_protein_from_het),
            'components': nx.number_connected_components(G_protein_from_het),
            'avg_degree': np.mean([d for n, d in G_protein_from_het.degree()]) if G_protein_from_het.number_of_nodes() > 0 else 0,
            'diameter': 0,
            'avg_path_length': 0
        }
    }
    
    # Calculate path-based metrics for connected components
    graphs = [
        ('Protein-only Network', G_protein),
        ('Heterogeneous Network (full)', G_het),
        ('Protein subgraph from Het.', G_protein_from_het)
    ]
    
    for name, graph in graphs:
        if nx.is_connected(graph):
            try:
                metrics_comparison[name]['diameter'] = nx.diameter(graph)
                metrics_comparison[name]['avg_path_length'] = nx.average_shortest_path_length(graph)
            except:
                pass
        else:
            # For disconnected graphs, analyze largest component
            if graph.number_of_nodes() > 0:
                largest_cc = max(nx.connected_components(graph), key=len)
                largest_subgraph = graph.subgraph(largest_cc)
                try:
                    metrics_comparison[name]['diameter'] = nx.diameter(largest_subgraph)
                    metrics_comparison[name]['avg_path_length'] = nx.average_shortest_path_length(largest_subgraph)
                except:
                    pass
    
    # Create comparison visualization
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    
    network_names = list(metrics_comparison.keys())
    colors = ['blue', 'red', 'green']
    
    metric_names = ['nodes', 'edges', 'density', 'avg_clustering', 'components', 'avg_degree', 'diameter', 'avg_path_length']
    
    for i, metric in enumerate(metric_names):
        ax = axes[i // 4, i % 4]
        values = [metrics_comparison[net][metric] for net in network_names]
        bars = ax.bar(network_names, values, color=colors)
        ax.set_title(f"{metric.replace('_', ' ').title()}")
        ax.tick_params(axis='x', rotation=45)
        
        # Add value labels on bars
        for bar, value in zip(bars, values):
            if value > 0:  # Only add label if value is meaningful
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(values)*0.01,
                       f'{value:.3f}' if isinstance(value, float) else str(value),
                       ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/network_approach_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save detailed comparison report
    with open(f"{output_path}/network_approach_comparison.txt", 'w') as f:
        f.write("NETWORK APPROACH COMPARISON\n")
        f.write("==========================\n\n")
        
        f.write("DETAILED METRICS COMPARISON:\n\n")
        
        for network_name, metrics in metrics_comparison.items():
            f.write(f"{network_name}:\n")
            for metric, value in metrics.items():
                f.write(f"  {metric}: {value:.6f if isinstance(value, float) else value}\n")
            f.write("\n")
        
        f.write("DECISION FRAMEWORK: WHEN TO USE WHICH APPROACH\n")
        f.write("=============================================\n\n")
        
        f.write("USE PROTEIN-ONLY NETWORK WHEN:\n")
        f.write("✓ Focus is on sequence evolution patterns\n")
        f.write("✓ Computational resources are limited\n")
        f.write("✓ Simple mutation pathway analysis\n")
        f.write("✓ Quick prototyping or exploratory analysis\n")
        f.write("✓ Dataset has limited non-protein annotations\n\n")
        
        f.write("USE HETEROGENEOUS NETWORK WHEN:\n")
        f.write("✓ Phenotype prediction is the main goal\n")
        f.write("✓ Rich metadata is available (organisms, ontologies, etc.)\n")
        f.write("✓ Complex machine learning models are planned\n")
        f.write("✓ Interpretability across biological layers is needed\n")
        f.write("✓ Transfer learning to other protein families\n\n")
        
        f.write("ALGORITHM RECOMMENDATIONS:\n")
        f.write("========================\n\n")
        
        f.write("FOR PROTEIN-ONLY NETWORKS:\n")
        f.write("- Node classification: GCN, GraphSAGE, GAT\n")
        f.write("- Community detection: Louvain, Leiden\n")
        f.write("- Centrality analysis: Betweenness, Eigenvector\n")
        f.write("- Motif analysis: Triangles, k-cores\n\n")
        
        f.write("FOR HETEROGENEOUS NETWORKS:\n")
        f.write("- Node classification: HGT, R-GCN, HAN\n")
        f.write("- Link prediction: HetSANN, R-GCN\n")
        f.write("- Graph-level tasks: HGP-SL, HGT\n")
        f.write("- Meta-path analysis: PathSim, HERec\n\n")
        
        # Add performance recommendations based on metrics
        protein_only_nodes = metrics_comparison['Protein-only Network']['nodes']
        het_nodes = metrics_comparison['Heterogeneous Network (full)']['nodes']
        
        f.write("PERFORMANCE CONSIDERATIONS:\n")
        f.write("=========================\n\n")
        f.write(f"Protein-only network: {protein_only_nodes} nodes\n")
        f.write(f"Heterogeneous network: {het_nodes} nodes\n")
        f.write(f"Size increase factor: {het_nodes/protein_only_nodes if protein_only_nodes > 0 else 0:.1f}x\n\n")
        
        if het_nodes / protein_only_nodes > 10:
            f.write("⚠️  RECOMMENDATION: Consider protein-only for initial analysis due to large size difference\n")
        elif het_nodes / protein_only_nodes > 3:
            f.write("📊 RECOMMENDATION: Use heterogeneous for final models, protein-only for prototyping\n")
        else:
            f.write("✅ RECOMMENDATION: Heterogeneous network is manageable, use for all analyses\n")
    
    LOGGER.info("Network approach comparison completed")
    return metrics_comparison


def main():
    """Main function to run both protein-only and heterogeneous network analyses"""
    LOGGER = setup_logging()
    LOGGER.info("Starting comprehensive network analysis (both approaches)")
    
    # Setup
    output_path = create_output_directory()
    eedb = setup_database_connection()
    
    # Analyze database schema first
    LOGGER.info("Analyzing database schema...")
    schema_info = analyze_database_schema(eedb)
    
    # Get TEM family data
    LOGGER.info("Fetching TEM family data...")
    data_tem_ids = get_tem_family_data(eedb)
    
    # ==================== PROTEIN-ONLY NETWORK ANALYSIS ====================
    LOGGER.info("="*60)
    LOGGER.info("PROTEIN-ONLY NETWORK ANALYSIS")
    LOGGER.info("="*60)
    
    LOGGER.info("Fetching protein network data...")
    results = get_protein_network_data(eedb)
    
    LOGGER.info("Creating protein-only NetworkX graph...")
    G_protein = create_networkx_graph(results)
    
    # Run protein-only analyses
    LOGGER.info("1. Analyzing node degrees (protein-only)...")
    degree_stats = analyze_node_degrees(G_protein, output_path)
    
    LOGGER.info("2. Analyzing betweenness centrality (protein-only)...")
    centrality_stats = analyze_betweenness_centrality(G_protein, output_path, data_tem_ids)
    
    LOGGER.info("3. Analyzing graph motifs (protein-only)...")
    motif_stats = analyze_graph_motifs(G_protein, output_path)
    
    LOGGER.info("4. Analyzing long distance interference (protein-only)...")
    distance_stats = analyze_long_distance_interference(G_protein, output_path)
    
    LOGGER.info("5. Analyzing mutation lengths by phenotype (protein-only)...")
    phenotype_stats = analyze_mutation_lengths_by_phenotype(eedb, results, output_path)
    
    # ==================== HETEROGENEOUS NETWORK ANALYSIS ====================
    LOGGER.info("="*60)
    LOGGER.info("HETEROGENEOUS NETWORK ANALYSIS")
    LOGGER.info("="*60)
    
    LOGGER.info("Fetching heterogeneous network data...")
    het_relationships = get_heterogeneous_network_data(eedb)
    
    if het_relationships:
        LOGGER.info("Creating heterogeneous NetworkX graph...")
        G_het = create_heterogeneous_networkx_graph(het_relationships)
        
        LOGGER.info("6. Analyzing heterogeneous network structure...")
        het_stats = analyze_heterogeneous_network_structure(G_het, output_path)
        
        # ==================== COMPARISON ANALYSIS ====================
        LOGGER.info("="*60)
        LOGGER.info("COMPARING BOTH APPROACHES")
        LOGGER.info("="*60)
        
        LOGGER.info("7. Comparing network approaches...")
        comparison_stats = compare_network_approaches(G_protein, G_het, output_path, data_tem_ids)
    else:
        LOGGER.warning("No heterogeneous relationships found. Skipping heterogeneous analysis.")
        G_het = None
        het_stats = {}
        comparison_stats = {}
    
    # ==================== COMPREHENSIVE SUMMARY ====================
    LOGGER.info("Creating comprehensive summary...")
    
    with open(f"{output_path}/comprehensive_network_analysis_summary.txt", 'w') as f:
        f.write("COMPREHENSIVE NETWORK ANALYSIS SUMMARY\n")
        f.write("=====================================\n\n")
        f.write(f"Analysis completed on: {pd.Timestamp.now()}\n")
        f.write(f"Output directory: {output_path}\n\n")
        
        f.write("DATABASE SCHEMA OVERVIEW:\n")
        f.write(f"Available node types: {len(schema_info['node_labels'])}\n")
        for label in schema_info['node_labels']:
            count = schema_info['node_counts'].get(label, 0)
            f.write(f"  - {label}: {count} nodes\n")
        f.write(f"Available relationship types: {len(schema_info['relationship_types'])}\n")
        for rel_type in schema_info['relationship_types']:
            f.write(f"  - {rel_type}\n")
        f.write("\n")
        
        f.write("PROTEIN-ONLY NETWORK OVERVIEW:\n")
        f.write(f"  Nodes: {G_protein.number_of_nodes()}\n")
        f.write(f"  Edges: {G_protein.number_of_edges()}\n")
        f.write(f"  Density: {nx.density(G_protein):.6f}\n")
        f.write(f"  Connected: {nx.is_connected(G_protein)}\n\n")
        
        if G_het:
            f.write("HETEROGENEOUS NETWORK OVERVIEW:\n")
            f.write(f"  Total nodes: {G_het.number_of_nodes()}\n")
            f.write(f"  Total edges: {G_het.number_of_edges()}\n")
            f.write(f"  Node types: {len(het_stats.get('node_type_counts', {}))}\n")
            f.write(f"  Relationship types: {len(het_stats.get('edge_type_counts', {}))}\n")
            f.write(f"  Density: {nx.density(G_het):.6f}\n\n")
        
        f.write("KEY FINDINGS (PROTEIN-ONLY):\n")
        f.write(f"  Mean degree: {degree_stats['Mean']:.3f}\n")
        f.write(f"  Network diameter: {distance_stats['diameter']}\n")
        f.write(f"  Average path length: {distance_stats['avg_path_length']:.3f}\n")
        f.write(f"  Small-world sigma: {distance_stats['small_world_sigma']:.3f}\n")
        f.write(f"  Number of triangles: {motif_stats['Triangles']}\n")
        f.write(f"  Number of analyzed phenotypes: {len(phenotype_stats)}\n\n")
        
        f.write("RECOMMENDATIONS FOR GRAPH ML ALGORITHM SELECTION:\n")
        f.write("================================================\n\n")
        
        protein_nodes = G_protein.number_of_nodes()
        het_nodes = G_het.number_of_nodes() if G_het else 0
        
        if het_nodes == 0:
            f.write("➤ Use protein-only approach (heterogeneous data not available)\n")
            f.write("  Recommended algorithms: GCN, GraphSAGE, GAT\n")
        elif het_nodes / protein_nodes < 2:
            f.write("➤ Use heterogeneous approach (manageable size increase)\n")
            f.write("  Recommended algorithms: HGT, R-GCN, HAN\n")
        elif het_nodes / protein_nodes < 5:
            f.write("➤ Use heterogeneous for final models, protein-only for prototyping\n")
            f.write("  Recommended algorithms: Start with GraphSAGE, scale to HGT\n")
        else:
            f.write("➤ Start with protein-only approach due to computational complexity\n")
            f.write("  Recommended algorithms: GraphSAGE, FastGCN for scalability\n")
        
        f.write("\nFILES GENERATED:\n")
        f.write("  PROTEIN-ONLY ANALYSIS:\n")
        f.write("    - degree_analysis.png\n")
        f.write("    - betweenness_centrality_analysis.png\n")
        f.write("    - motif_analysis.png\n")
        f.write("    - long_distance_interference_analysis.png\n")
        f.write("    - mutation_length_phenotype_analysis.png\n")
        if G_het:
            f.write("  HETEROGENEOUS ANALYSIS:\n")
            f.write("    - heterogeneous_network_analysis.png\n")
            f.write("  COMPARISON:\n")
            f.write("    - network_approach_comparison.png\n")
        f.write("  STATISTICS:\n")
        f.write("    - Various .txt statistics files\n")
    
    LOGGER.info(f"Analysis completed successfully. Results saved to: {output_path}")
    print(f"\nComprehensive network analysis completed!")
    print(f"Results saved to: {output_path}")
    print(f"Protein-only network: {G_protein.number_of_nodes()} nodes, {G_protein.number_of_edges()} edges")
    if G_het:
        print(f"Heterogeneous network: {G_het.number_of_nodes()} nodes, {G_het.number_of_edges()} edges")
        print(f"Recommendation: {'Heterogeneous' if het_nodes/protein_nodes < 3 else 'Start with protein-only'} approach")


if __name__ == "__main__":
    main()


    # nohup python scr/code/022_Plots_Network_Analysis_Staab.py > 022_Plots_Network_Analysis_Staab.log 2>&1 &
