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


def get_protein_network_data(eedb, name_of_standard_numbering_tool="standard_numbering_pairwise_flagged_proteins", exclude_synthetic=True, allow_multi_mutations=False, include_all_lengths=False):
    """Query Neo4j to get protein mutation network data
    
    Parameters:
    -----------
    eedb : Pyeed
        Database connection object
    name_of_standard_numbering_tool : str
        Name of the standard numbering tool to use
    exclude_synthetic : bool
        If True, exclude synthetic organisms (taxonomy_id = 32630)
        If False, include all organisms
    allow_multi_mutations : bool
        If True, include proteins connected by multiple mutations
        If False, only include single mutations (default)
    include_all_lengths : bool
        If True, include proteins of all sequence lengths
        If False, only include proteins with length 286 (default)
    """
    LOGGER = logging.getLogger(__name__)
    
    # Base query parts
    synthetic_filter = "AND o1.taxonomy_id <> 32630 AND o2.taxonomy_id <> 32630" if exclude_synthetic else ""
    
    # Build conditional filters
    length_filter = "AND p1.seq_length = 286 AND p2.seq_length = 286" if not include_all_lengths else ""
    mutation_filter = "AND size(r.from_positions) = 1" if not allow_multi_mutations else ""
    
    query = f"""
    MATCH (p1:Protein)-[r:MUTATION]->(p2:Protein)
    WHERE 1=1 {length_filter} {mutation_filter}
    WITH p1, p2, r
    MATCH (o1:Organism)-[:ORIGINATES_FROM]-(p1)
    MATCH (o2:Organism)-[:ORIGINATES_FROM]-(p2)
    WHERE 1=1 {synthetic_filter}
    WITH p1, p2, r, o1, o2
    MATCH (s:StandardNumbering {{name: $name_of_standard_numbering_tool}})-[s1:HAS_STANDARD_NUMBERING]-(p1)
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
           o2.name AS target_organism,
           o1.taxonomy_id AS source_taxonomy_id,
           o2.taxonomy_id AS target_taxonomy_id
    """

    results = eedb.db.execute_read(query, parameters={"name_of_standard_numbering_tool": name_of_standard_numbering_tool})
    
    # Build descriptive log message
    organism_type = "excluding synthetic" if exclude_synthetic else "including synthetic"
    mutation_type = "multi-mutations" if allow_multi_mutations else "single mutations only"
    length_type = "all lengths" if include_all_lengths else "length 286 only"
    
    LOGGER.info(f"Number of mutation relationships ({organism_type}, {mutation_type}, {length_type}): {len(results)}")
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


def get_heterogeneous_network_data(eedb, exclude_synthetic=True):
    """Query Neo4j to get heterogeneous network data including multiple node types
    
    Parameters:
    -----------
    eedb : Pyeed
        Database connection object
    exclude_synthetic : bool
        If True, exclude synthetic organisms (taxonomy_id = 32630)
        If False, include all organisms
    """
    LOGGER = logging.getLogger(__name__)
    
    # Base filter for synthetic organisms
    synthetic_filter = "AND o.taxonomy_id <> 32630" if exclude_synthetic else ""
    synthetic_filter_protein = "AND o1.taxonomy_id <> 32630 AND o2.taxonomy_id <> 32630" if exclude_synthetic else ""
    
    # Query for Protein-Protein MUTATION relationships
    protein_mutation_query = f"""
    MATCH (p1:Protein)-[r:MUTATION]->(p2:Protein)
    WHERE p1.seq_length = 286 AND p2.seq_length = 286
    {"WITH p1, p2, r MATCH (o1:Organism)-[:ORIGINATES_FROM]-(p1) MATCH (o2:Organism)-[:ORIGINATES_FROM]-(p2) WHERE 1=1" + synthetic_filter_protein + " WITH p1, p2, r" if exclude_synthetic else ""}
    RETURN 
        'Protein' as source_type, p1.accession_id AS source_id,
        'Protein' as target_type, p2.accession_id AS target_id,
        'MUTATION' as relationship_type,
        r.from_positions as mutation_positions
    """
    
    # Query for Protein-Organism relationships
    protein_organism_query = f"""
    MATCH (p:Protein)-[r:ORIGINATES_FROM]->(o:Organism)
    WHERE p.seq_length = 286 {synthetic_filter}
    RETURN 
        'Protein' as source_type, p.accession_id AS source_id,
        'Organism' as target_type, o.taxonomy_id AS target_id,
        'ORIGINATES_FROM' as relationship_type,
        o.name as organism_name
    """
    
    # Query for DNA-Protein relationships
    dna_protein_query = f"""
    MATCH (d:DNA)-[r:ENCODES]->(p:Protein)
    WHERE p.seq_length = 286
    {"WITH d, r, p MATCH (o:Organism)-[:ORIGINATES_FROM]-(p) WHERE 1=1" + synthetic_filter + " WITH d, r, p" if exclude_synthetic else ""}
    RETURN 
        'DNA' as source_type, d.accession_id AS source_id,
        'Protein' as target_type, p.accession_id AS target_id,
        'ENCODES' as relationship_type,
        d.gc_content as gc_content
    """
    
    # Query for Protein-OntologyObject relationships (resistance mechanisms)
    protein_ontology_query = f"""
    MATCH (p:Protein)-[r]-(o:OntologyObject)
    WHERE p.seq_length = 286
    {"WITH p, r, o MATCH (org:Organism)-[:ORIGINATES_FROM]-(p) WHERE 1=1" + synthetic_filter.replace("o.taxonomy_id", "org.taxonomy_id") + " WITH p, r, o" if exclude_synthetic else ""}
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
    
    organism_type = "excluding synthetic" if exclude_synthetic else "including synthetic"
    
    for query_name, query in queries:
        try:
            results = eedb.db.execute_read(query)
            LOGGER.info(f"{query_name} ({organism_type}): {len(results)} relationships")
            all_relationships.extend(results)
        except Exception as e:
            LOGGER.warning(f"Query {query_name} failed: {e}")
    
    LOGGER.info(f"Total heterogeneous relationships collected ({organism_type}): {len(all_relationships)}")
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


def save_network_info(G, output_path, network_type, exclude_synthetic=True):
    """Save basic network information to file"""
    organism_suffix = "_excl_synthetic" if exclude_synthetic else "_incl_synthetic"
    filename = f"{output_path}/network_info_{network_type}{organism_suffix}.txt"
    
    with open(filename, 'w') as f:
        f.write(f"NETWORK INFORMATION: {network_type.upper()}\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Organism inclusion: {'Excluding synthetic' if exclude_synthetic else 'Including synthetic'}\n")
        f.write(f"Number of nodes: {G.number_of_nodes()}\n")
        f.write(f"Number of edges: {G.number_of_edges()}\n")
        f.write(f"Network density: {nx.density(G):.6f}\n")
        f.write(f"Is connected: {nx.is_connected(G)}\n")
        f.write(f"Number of connected components: {nx.number_connected_components(G)}\n")
        
        if G.number_of_nodes() > 0:
            largest_cc_size = len(max(nx.connected_components(G), key=len))
            f.write(f"Largest connected component size: {largest_cc_size}\n")
            f.write(f"Largest component fraction: {largest_cc_size / G.number_of_nodes():.3f}\n")
