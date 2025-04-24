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
password = os.getenv("NEO4J_NIKLAS_TEM_CLEAN")
if password is None:
    raise ValueError("KEY is not set in the .env file.")

uri = "bolt://129.69.129.130:2123"
user = "neo4j"
eedb = Pyeed(uri, user=user, password=password)
eedb.db.initialize_db_constraints(user, password)

name_of_standard_numbering_tool = (
    "standard_numbering_pairwise_circular_mutations_to_blaTEM1a"
)
et = EmbeddingTool()
sn = StandardNumberingTool(name=name_of_standard_numbering_tool)
md = MutationDetection()

blaTEM1a_database_id = "CAD09800.1"
max_number_of_mutations = 10


if __name__ == "__main__":
    # Read and prepare data
    path = "/home/nab/Niklas/TEM-lactamase/combined_protein_data.csv"
    df = pd.read_csv(path)

    names = [
        "TEM beta-lactamase",
        "AER beta-lactamase",
        "DES beta-lactamase",
        "CKO beta-lactamase",
        "BKC Beta-lactamase",
        "BIC Beta-lactamase",
        "HERA beta-lactamase",
        "GES beta-lactamase",
    ]

    df_tem = df[df["family_name"].isin(names)]
    ids_in_circle_series = df_tem["ids_in_circle"].apply(
        lambda x: ast.literal_eval(x) if x != "[]" else []
    )
    ids_in_circle = []
    for circle_list in ids_in_circle_series:
        ids_in_circle.extend(circle_list)
    ids_in_circle = list(dict.fromkeys(ids_in_circle))

    LOGGER.info(f"Number of IDs in circle: {len(ids_in_circle)}")

    # Query to find proteins with exactly one mutation between them
    query = """
    MATCH (p:Protein)-[r:MUTATION]-(q:Protein)
    WHERE p.accession_id IN $ids_in_circle AND q.accession_id IN $ids_in_circle AND size(r.from_positions) = 1
    RETURN p.accession_id AS source, q.accession_id AS target, r.from_positions[0] AS from_pos, r.to_positions[0] AS to_pos
    LIMIT 4000
    """

    results = eedb.db.execute_read(query, parameters={"ids_in_circle": ids_in_circle})
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
    plt.figure(figsize=(15, 15))
    pos = nx.spring_layout(G, k=0.3, iterations=50)

    # make a betweenness centrality plot
    betweenness_centrality = nx.betweenness_centrality(G)
    # plot the betweenness centrality
    plt.figure(figsize=(10, 8))
    plt.bar(betweenness_centrality.keys(), betweenness_centrality.values())
    plt.title("Betweenness Centrality")
    plt.xlabel("Node")
    plt.ylabel("Betweenness Centrality")
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
