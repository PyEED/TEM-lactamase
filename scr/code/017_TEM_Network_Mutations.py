import ast
import logging
import os
from collections import Counter

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.embedding_analysis import EmbeddingTool
from pyeed.analysis.mutation_detection import MutationDetection
from pyeed.analysis.standard_numbering import StandardNumberingTool

path_to_data_blast = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-01-19_12-37-48"

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)


load_dotenv()
password = os.getenv("NEO4J_NIKLAS_TEM_THREE")
if password is None:
    raise ValueError("KEY is not set in the .env file.")

uri = "bolt://129.69.129.130:2137"
user = "neo4j"
eedb = Pyeed(uri, user=user, password=password)
# eedb.db.initialize_db_constraints(user, password)

name_of_standard_numbering_tool = (
    "standard_numbering_pairwise_flagged_proteins"
)
et = EmbeddingTool()
sn = StandardNumberingTool(name=name_of_standard_numbering_tool)
md = MutationDetection()

blaTEM1a_database_id = "CAD09800.1"
max_number_of_mutations = 10


if __name__ == "__main__":
    # Query to find proteins with exactly one mutation between them
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
    MATCH (s:StandardNumbering {name: 'standard_numbering_pairwise_flagged_proteins'})-[s1:HAS_STANDARD_NUMBERING]-(p1)
    MATCH (s)-[s2:HAS_STANDARD_NUMBERING]-(p2)
    WITH p1, p2, r, s1, s2, o1, o2, rand() AS random_value, id(r) AS id_r
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

    results = eedb.db.execute_read(query)
    LOGGER.info(f"Number of results: {len(results)}")

    # Create an undirected graph
    G = nx.Graph()

    # Add nodes and edges
    for record in results:
        source = record["source"]
        target = record["target"]
        from_pos = record["from_pos"]
        to_pos = record["to_pos"]

        # Add nodes if they don't exist
        if source not in G:
            G.add_node(source)
        if target not in G:
            G.add_node(target)

        # Add edge with mutation information
        G.add_edge(source, target, from_pos=from_pos, to_pos=to_pos)

    # make a degree distribution plot
    degree_sequence = [d for n, d in G.degree()]
    degree_counts = Counter(degree_sequence)
    degrees = sorted(degree_counts.keys())
    counts = [degree_counts[d] for d in degrees]

    plt.figure(figsize=(10, 8))
    plt.loglog(degrees, counts, "bo-", linewidth=2, markersize=8)
    plt.grid(True, which="both", ls="-", alpha=0.6)
    plt.title("Degree Distribution (Log-Log Scale)")
    plt.xlabel("Degree (log)")
    plt.ylabel("Frequency (log)")
    plt.tight_layout()
    plt.savefig("degree_distribution.png", dpi=300)
    plt.close()

    # Create visualization
    plt.figure(figsize=(30, 30))
    pos = nx.spring_layout(G, k=2.0, iterations=200, scale=2.0)

    # make a betweenness centrality plot
    betweenness_centrality = nx.betweenness_centrality(G)
    # Get the top 10 nodes by betweenness centrality
    top_nodes = sorted(
        betweenness_centrality.items(), key=lambda x: x[1], reverse=True
    )[:10]
    top_nodes_dict = {node: value for node, value in top_nodes}

    # plot the betweenness centrality for top 10 nodes
    plt.figure(figsize=(10, 8))
    plt.bar(top_nodes_dict.keys(), top_nodes_dict.values())
    plt.title("Betweenness Centrality (Top 10 Nodes)")
    plt.xlabel("Node")
    plt.ylabel("Betweenness Centrality")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig("betweenness_centrality.png", dpi=300)
    plt.close()

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=100, node_color="lightblue")

    # Draw edges without labels
    nx.draw_networkx_edges(G, pos, arrows=False)

    plt.title("Protein Mutation Network (Single Mutations Only)")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig("protein_mutation_network.png", dpi=300)
    plt.close()

    LOGGER.info(
        f"Network created with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges"
    )
    LOGGER.info("Network visualization saved as 'protein_mutation_network.png'")

    # nohup python scr/code/017_TEM_Network_Mutations.py > 017_TEM_Network_Mutations.log 2>&1 &
