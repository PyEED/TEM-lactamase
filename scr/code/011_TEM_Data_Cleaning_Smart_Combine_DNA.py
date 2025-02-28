# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import json
import logging
import os

from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.embedding_analysis import EmbeddingTool

# ------------------------------------- SETUP -------------------------------------


path_to_data_blast_dna = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-01-19_12-37-48"
path_to_data_blast_protein = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data/combined_data_blast_5000_tem_209"


load_dotenv()
password = os.getenv("NEO4J_NIKLAS_TEM_CLEAN")
if password is None:
    raise ValueError("KEY is not set in the .env file.")


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)


uri = "bolt://129.69.129.130:2123"
user = "neo4j"
eedb = Pyeed(uri, user=user, password=password)
eedb.db.initialize_db_constraints(user, password)

et = EmbeddingTool()

# ------------------------------------- FUNCTIONS -------------------------------------


def get_all_realtiontship_for_node(node_id, node_type):
    """
    Get all relationships for a node. So they can be recreated for the other node that is kept.
    """
    query = f"""
    MATCH (n:{node_type} {{accession_id: "{node_id}"}})-[r]-(e) 
    RETURN type(r) as type, properties(r) as properties, 
           startNode(r) as start_node, id(startNode(r)) as start_node_id,
           endNode(r) as end_node, id(endNode(r)) as end_node_id,
           labels(startNode(r)) as start_labels, head(labels(startNode(r))) as start_node_type,
           labels(endNode(r)) as end_labels, head(labels(endNode(r))) as end_node_type,
           properties(startNode(r)) as start_properties, properties(endNode(r)) as end_properties
    """
    return eedb.db.execute_read(query)


# ---------------------- Newly Outsourced Functions ---------------------- #


def remove_dna_node(dna_node_to_remove, dna_node_kept, eedb, logger):
    """
    Remove a DNA node and take over all the relationships of the removed one to the DNA node which is kept.
    """
    # Retrieve all relationships of the removed DNA node.
    relationships = get_all_realtiontship_for_node(
        dna_node_to_remove["d.accession_id"], "DNA"
    )

    # Remove the removed DNA node from the database.
    query_remove = f"""
        MATCH (d:DNA) WHERE d.accession_id = "{dna_node_to_remove['d.accession_id']}" DETACH DELETE d
    """
    eedb.db.execute_write(query_remove)
    logger.info(f"Removing {dna_node_to_remove['d.accession_id']} from the database")

    # Recreate the removed DNA node's relationships, reassigning to current DNA node.
    for relationship in relationships:
        start_node_label = relationship["start_node_type"]
        end_node_label = relationship["end_node_type"]

        # Determine which side in the relationship is the deleted (identical) protein.
        if (
            relationship["start_node"].get("accession_id")
            == dna_node_to_remove["d.accession_id"]
        ):
            # Removed DNA node is the start node.
            new_start_info = {
                "match_by": "accession",
                "value": dna_node_kept["d.accession_id"],
                "label": start_node_label,
            }
            if "accession_id" in relationship["end_node"]:
                new_end_info = {
                    "match_by": "accession",
                    "value": relationship["end_node"]["accession_id"],
                    "label": end_node_label,
                }
            else:
                new_end_info = {
                    "match_by": "id",
                    "value": relationship["end_node_id"],
                    "label": end_node_label,
                }
        elif (
            relationship["end_node"].get("accession_id")
            == dna_node_to_remove["d.accession_id"]
        ):
            # Removed DNA node is the end node.
            new_end_info = {
                "match_by": "accession",
                "value": dna_node_kept["d.accession_id"],
                "label": end_node_label,
            }
            if "accession_id" in relationship["start_node"]:
                new_start_info = {
                    "match_by": "accession",
                    "value": relationship["start_node"]["accession_id"],
                    "label": start_node_label,
                }
            else:
                new_start_info = {
                    "match_by": "id",
                    "value": relationship["start_node_id"],
                    "label": start_node_label,
                }
        else:
            logger.error(
                f"Relationship does not involve identical protein: {relationship}"
            )
            continue

        # Build the MATCH clauses dynamically based on whether to use property or internal id.
        if new_start_info["match_by"] == "accession":
            start_match_clause = f'MATCH (start:{new_start_info["label"]}) WHERE start.accession_id = "{new_start_info["value"]}"'
        else:
            start_match_clause = f'MATCH (start:{new_start_info["label"]}) WHERE id(start) = {new_start_info["value"]}'

        if new_end_info["match_by"] == "accession":
            end_match_clause = f'MATCH (end:{new_end_info["label"]}) WHERE end.accession_id = "{new_end_info["value"]}"'
        else:
            end_match_clause = f'MATCH (end:{new_end_info["label"]}) WHERE id(end) = {new_end_info["value"]}'

        # New code: convert any relationship properties to a Cypher map literal.
        rel_properties = relationship.get("properties", {})
        if rel_properties:
            # Construct a string like 'key1: "value1", key2: 5'
            props_str = ", ".join(
                f"{key}: {json.dumps(value)}" for key, value in rel_properties.items()
            )
            properties_literal = f"{{{props_str}}}"
        else:
            properties_literal = ""

        query_recreate_relationship = f"""
            {start_match_clause}
            {end_match_clause}
            CREATE (start)-[:{relationship['type']}{properties_literal}]->(end)
        """
        print(
            f"Recreating relationship of type {relationship['type']} from {new_start_info['value']} to {new_end_info['value']}"
        )
        logger.info(query_recreate_relationship)
        eedb.db.execute_write(query_recreate_relationship)


if __name__ == "__main__":
    # Here we are interested in starting a data cleaning.
    # We assume the proteins are already cleaned and combined. But there can be multiple DNA sequences for the same protein. Those might be identical. And we want to combine them.
    # As with the protein cleaning, we want to combine the double DNA sequences in the attribute 'IdenticalIds'.
    # We then also want to remove the duplicate DNA nodes. But keep for the DNA node that is kept all the realtionship of the removed ones.
    # Instead of searching with the vector index we just go through the list of proteins and check what DNA sequences conntect to them. And then one by one check wether they are identical in the relevant area aka the start and end of the DNA sequence.

    query_protein_ids = """
        MATCH (p:Protein) RETURN p.accession_id
    """
    protein_ids = eedb.db.execute_read(query_protein_ids)
    protein_ids = [protein_id["p.accession_id"] for protein_id in protein_ids]
    print(f"Number of proteins: {len(protein_ids)}")

    # Process each protein.
    for index in range(0, len(protein_ids)):
        current_protein_id = protein_ids[index]

        # Get all DNA nodes connected to the protein. The relationship goes from DNA to Protein and is called 'ENCODES', it goes from DNA to Protein.
        query_dna_nodes = f"""
            MATCH (d:DNA)-[r:ENCODES]->(p:Protein) WHERE p.accession_id = "{current_protein_id}" RETURN d.accession_id, d.sequence, r.start, r.end
        """
        dna_nodes = eedb.db.execute_read(query_dna_nodes)
        dna_nodes_to_remove = []

        for i in range(0, len(dna_nodes)):
            dna_node = dna_nodes[i]
            for j in range(i + 1, len(dna_nodes)):
                other_dna_node = dna_nodes[j]

                # Check if the DNA sequences are identical in the relevant area.
                if (
                    dna_node["r.start"] == other_dna_node["r.start"]
                    and dna_node["r.end"] == other_dna_node["r.end"]
                ):
                    LOGGER.info(
                        f"Identical DNA sequences: {dna_node['d.accession_id']} and {other_dna_node['d.accession_id']}"
                    )

                    # add the other DNA node to the list of DNA nodes to remove.
                    dna_nodes_to_remove.append((dna_node, other_dna_node))

        # remove the DNA nodes to remove.
        for dna_node_kept, dna_node_to_remove in dna_nodes_to_remove:
            remove_dna_node(dna_node_to_remove, dna_node_kept, eedb, LOGGER)
