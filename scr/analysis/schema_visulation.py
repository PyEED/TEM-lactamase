
from py2neo import Graph
import networkx as nx
import matplotlib.pyplot as plt

# Connect to Neo4j (adjust connection details as necessary)
graph = Graph("bolt://localhost:7687", auth=("neo4j", "12345678"))

# Run the Cypher query to get schema information
query = """
MATCH (n)-[r]->(m)
RETURN DISTINCT labels(n) AS source_label, type(r) AS relationship, labels(m) AS target_label
"""
results = graph.run(query).data()

# Initialize a directed graph with NetworkX
G = nx.DiGraph()

# Process the results and add nodes and edges to the graph
for record in results:
    source_label = record['source_label'][0]  # Get the first label from the list
    target_label = record['target_label'][0]  # Get the first label from the list
    relationship = record['relationship']
    
    # Add nodes and edges to the NetworkX graph
    G.add_node(source_label, label=source_label)
    G.add_node(target_label, label=target_label)
    G.add_edge(source_label, target_label, relationship=relationship)

# Draw the graph
plt.figure(figsize=(10, 10))
# which layout: spring_layout, shell_layout, random_layout, circular_layout, kamada_kawai_layout, planar_layout, spectral_layout
pos = nx.planar_layout(G)
edge_labels = {(u, v): d['relationship'] for u, v, d in G.edges(data=True)}
# nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color="red", rotate=True)
# I want big arrow heads
nx.draw_networkx_edges(G, pos, connectionstyle="arc3,rad=0.1", edge_color="black", arrowsize=10)
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color="gray", label_pos=0.4, font_size=10)
nx.draw(G, pos, with_labels=True, node_color="lightgray", node_size=5000, font_size=10, font_weight="bold", edgelist=[])

# Save and show the plot
plt.savefig("neo4j_schema_visualization.png", format="PNG")
plt.show()