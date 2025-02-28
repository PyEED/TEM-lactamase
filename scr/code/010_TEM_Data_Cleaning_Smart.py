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
password = os.getenv("NEO4J_NIKLAS_TEM_HARRY")
if password is None:
    raise ValueError("KEY is not set in the .env file.")


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)


uri = "bolt://129.69.129.130:4123"
user = "neo4j"
eedb = Pyeed(uri, user=user, password=password)
eedb.db.initialize_db_constraints(user, password)

et = EmbeddingTool()

# ------------------------------------- FUNCTIONS -------------------------------------


def label_node(node_type, node_id, node_label, label_type):
    """
    Label a node in the Neo4j database.
    The node type is the type of the node to be labeled. Could be Protein, DNA
    The node id is the id of the node to be labeled. Most likely a accession_id
    The node label is the label of the node to be labeled. Could be DNA_Blast, Protein_Blast, Protein_Bldb
    The label_type is the type of label to be added. Could be "Source", "TEM_Type", "Resistance_Mechanism", "Species"

    The label is an attribute of the node.
    """

    query = f"""
    MATCH (n:{node_type} {{accession_id: "{node_id}"}})
    SET n.{label_type} = "{node_label}"
    """

    eedb.db.execute_write(query)


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


def remove_and_reassign(current_protein_id, identical_protein_id, eedb, logger):
    """
    Remove the identical protein node and reassign all its relationships
    to the current protein node.

    Parameters:
        current_protein_id (str): The protein ID of the node to keep.
        identical_protein_id (str): The duplicate protein ID to be removed.
        eedb (Pyeed): The database connection object.
        logger (logging.Logger): Logger to record execution details.
    """
    # Retrieve all relationships of the identical protein node.
    relationships = get_all_realtiontship_for_node(identical_protein_id, "Protein")

    # Remove the identical protein node from the database.
    query_remove = f"""
        MATCH (p:Protein) WHERE p.accession_id = "{identical_protein_id}" DETACH DELETE p
    """
    eedb.db.execute_write(query_remove)
    logger.info(f"Removing {identical_protein_id} from the database")

    # Recreate the removed protein's relationships, reassigning to current protein node.
    for relationship in relationships:
        start_node_label = relationship["start_node_type"]
        end_node_label = relationship["end_node_type"]

        # Determine which side in the relationship is the deleted (identical) protein.
        if relationship["start_node"].get("accession_id") == identical_protein_id:
            # Removed protein is the start node.
            new_start_info = {
                "match_by": "accession",
                "value": current_protein_id,
                "label": "Protein",
            }
            # For the partner (end node), use accession_id if available; otherwise use the internal node id.
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
        elif relationship["end_node"].get("accession_id") == identical_protein_id:
            # Removed protein is the end node.
            new_end_info = {
                "match_by": "accession",
                "value": current_protein_id,
                "label": "Protein",
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
            properties_literal = f" {{{props_str}}}"
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


def update_and_remove_identical(current_protein_id, identical_protein_id, eedb, logger):
    """
    Update the IdenticalIds attribute of the current protein node with the identical protein's ID,
    then remove the identical protein node and reassign its relationships.

    Parameters:
        current_protein_id (str): The protein ID of the node to keep.
        identical_protein_id (str): The duplicate protein ID to be removed.
        eedb (Pyeed): The database connection object.
        logger (logging.Logger): Logger to record execution details.
    """
    # Retrieve the IdenticalIds from the kept protein node.
    query_get_kept_ids = f"""
        MATCH (p:Protein {{accession_id: "{current_protein_id}"}})
        RETURN COALESCE(p.IdenticalIds, []) AS IdenticalIds
    """
    result_kept = eedb.db.execute_read(query_get_kept_ids)
    kept_ids = result_kept[0].get("IdenticalIds", []) if result_kept else []

    # Retrieve the IdenticalIds from the removed protein node.
    query_get_removed_ids = f"""
        MATCH (p:Protein {{accession_id: "{identical_protein_id}"}})
        RETURN COALESCE(p.IdenticalIds, []) AS IdenticalIds
    """
    result_removed = eedb.db.execute_read(query_get_removed_ids)
    removed_ids = result_removed[0].get("IdenticalIds", []) if result_removed else []

    # Combine the kept and removed IdenticalIds and add the removed protein's ID.
    new_ids = set(kept_ids) | set(removed_ids)
    new_ids.add(identical_protein_id)
    new_ids_list = list(new_ids)

    query_write_identical_ids = f"""
        MATCH (p:Protein) WHERE p.accession_id = "{current_protein_id}" 
        SET p.IdenticalIds = {json.dumps(new_ids_list)}
    """
    eedb.db.execute_write(query_write_identical_ids)

    # Remove the identical protein and reassign its relationships.
    remove_and_reassign(current_protein_id, identical_protein_id, eedb, logger)


def handle_identical_protein(current_protein_id, identical_protein_id, eedb, logger):
    """
    Process two proteins that are identified as identical.
    The function checks the 'Source' attribute for both proteins to determine if they are DE_NOVO_BASED_ON_DNA.
    Depending on the combination:
        - If both proteins or only the identical protein is DE_NOVO_BASED_ON_DNA, the identical protein node is removed
          and its relationships reassigned.
        - If neither protein is DE_NOVO_BASED_ON_DNA, the identical protein's ID is added to the IdenticalIds of the
          current protein, and then the identical protein node is removed with relationship reassignment.

    Parameters:
        current_protein_id (str): The protein ID of the node to keep.
        identical_protein_id (str): The duplicate protein ID detected as identical.
        eedb (Pyeed): The database connection object.
        logger (logging.Logger): Logger to record execution details.
    """
    query_identical_protein = f"""
        MATCH (p:Protein) WHERE p.accession_id = "{identical_protein_id}" RETURN p.Source
    """
    identical_protein_source = eedb.db.execute_read(query_identical_protein)
    identical_protein_is_de_novo = (
        identical_protein_source[0]["p.Source"] == "DE_NOVO_BASED_ON_DNA"
    )

    query_current_protein = f"""
        MATCH (p:Protein) WHERE p.accession_id = "{current_protein_id}" RETURN p.Source
    """
    current_protein_source = eedb.db.execute_read(query_current_protein)
    current_protein_is_de_novo = (
        current_protein_source[0]["p.Source"] == "DE_NOVO_BASED_ON_DNA"
    )

    if identical_protein_is_de_novo and current_protein_is_de_novo:
        logger.info(
            f"Both {current_protein_id} and {identical_protein_id} are DE_NOVO_BASED_ON_DNA"
        )
        remove_and_reassign(current_protein_id, identical_protein_id, eedb, logger)
    elif identical_protein_is_de_novo and not current_protein_is_de_novo:
        logger.info(
            f"{current_protein_id} is not DE_NOVO_BASED_ON_DNA and {identical_protein_id} is DE_NOVO_BASED_ON_DNA"
        )
        remove_and_reassign(current_protein_id, identical_protein_id, eedb, logger)
    else:
        logger.info(
            f"{current_protein_id} and {identical_protein_id} are not DE_NOVO_BASED_ON_DNA"
        )
        update_and_remove_identical(
            current_protein_id, identical_protein_id, eedb, logger
        )


if __name__ == "__main__":
    # Here we are interested in starting a data cleaning.
    # First we want to find identical protein sequences.

    """
    et.drop_vector_index(index_name="vector_index_Protein_embedding", db=eedb.db)

    et.create_embedding_vector_index_neo4j(
        index_name="vector_index_Protein_embedding",
        db=eedb.db,
        similarity_function="cosine",
        m=512,
        ef_construction=3200,
        dimensions=960,
    )
    """

    query_protein_ids = """
        MATCH (p:Protein) RETURN p.accession_id
    """
    protein_ids = eedb.db.execute_read(query_protein_ids)
    protein_ids = [protein_id["p.accession_id"] for protein_id in protein_ids]
    print(f"Number of proteins: {len(protein_ids)}")

    # Process each protein.
    for index in range(0, len(protein_ids)):
        current_protein_id = protein_ids[index]
        results = et.find_nearest_neighbors_based_on_vector_index(
            index_name="vector_index_Protein_embedding",
            query_protein_id=current_protein_id,
            number_of_neighbors=10,
            db=eedb.db,
        )

        for result_list in results:
            for i in range(1, len(result_list)):
                if result_list[i] == 1.0:
                    # Outsource handling of the identical protein into its own function.
                    if current_protein_id != result_list[0]:
                        LOGGER.info(
                            f"Found one that is identical to {current_protein_id}"
                        )
                        handle_identical_protein(
                            current_protein_id, result_list[i], eedb, LOGGER
                        )
