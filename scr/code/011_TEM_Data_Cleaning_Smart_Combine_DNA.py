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
password = os.getenv("NEO4J_NIKLAS_TEM_THREE")
if password is None:
    raise ValueError("KEY is not set in the .env file.")

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)

uri = "bolt://129.69.129.130:2137"
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
    logger.info(
        f"Retrieving all relationships of {dna_node_to_remove['d.accession_id']}"
    )
    relationships = get_all_realtiontship_for_node(
        dna_node_to_remove["d.accession_id"], "DNA"
    )
    logger.info(f"Retrieved {len(relationships)} relationships")
    # Remove the removed DNA node from the database.
    query_remove = f"""
        MATCH (d:DNA) WHERE d.accession_id = "{dna_node_to_remove["d.accession_id"]}" DETACH DELETE d
    """
    eedb.db.execute_write(query_remove)
    logger.info(f"Removing {dna_node_to_remove['d.accession_id']} from the database")

    # Recreate the removed DNA node's relationships, reassigning to the DNA node which is kept.
    for relationship in relationships:
        start_node_label = relationship["start_node_type"]
        end_node_label = relationship["end_node_type"]

        # Determine which side in the relationship is the removed DNA node.
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

        # Build the MATCH clauses dynamically.
        if new_start_info["match_by"] == "accession":
            start_match_clause = f'MATCH (start:{new_start_info["label"]}) WHERE start.accession_id = "{new_start_info["value"]}"'
        else:
            start_match_clause = f"MATCH (start:{new_start_info['label']}) WHERE id(start) = {new_start_info['value']}"

        if new_end_info["match_by"] == "accession":
            end_match_clause = f'MATCH (end:{new_end_info["label"]}) WHERE end.accession_id = "{new_end_info["value"]}"'
        else:
            end_match_clause = f"MATCH (end:{new_end_info['label']}) WHERE id(end) = {new_end_info['value']}"

        # Convert relationship properties to a Cypher map literal.
        rel_properties = relationship.get("properties", {})
        if rel_properties:
            props_str = ", ".join(
                f"{key}: {json.dumps(value)}" for key, value in rel_properties.items()
            )
            properties_literal = f"{{{props_str}}}"
        else:
            properties_literal = ""

        query_recreate_relationship = f"""
            {start_match_clause}
            {end_match_clause}
            CREATE (start)-[:{relationship["type"]}{properties_literal}]->(end)
        """
        logger.info(query_recreate_relationship)
        print(
            f"Recreating relationship of type {relationship['type']} from {new_start_info['value']} to {new_end_info['value']}"
        )
        eedb.db.execute_write(query_recreate_relationship)


def update_and_remove_identical_dna(dna_node_kept, dna_node_to_remove, eedb, logger):
    """
    Update the IdenticalIds attribute of the kept DNA node (if it exists) with the removed DNA node's accession,
    then remove the duplicate DNA node and reassign its relationships.
    """
    # Retrieve the IdenticalIds from the kept DNA node.
    query_get_kept_ids = f"""
        MATCH (d:DNA {{accession_id: "{dna_node_kept["d.accession_id"]}"}})
        RETURN COALESCE(d.IdenticalIds, []) AS IdenticalIds
    """
    result_kept = eedb.db.execute_read(query_get_kept_ids)
    kept_ids = result_kept[0].get("IdenticalIds", []) if result_kept else []

    # Retrieve the IdenticalIds from the removed DNA node.
    query_get_removed_ids = f"""
        MATCH (d:DNA {{accession_id: "{dna_node_to_remove["d.accession_id"]}"}})
        RETURN COALESCE(d.IdenticalIds, []) AS IdenticalIds
    """
    result_removed = eedb.db.execute_read(query_get_removed_ids)
    removed_ids = result_removed[0].get("IdenticalIds", []) if result_removed else []

    # Combine kept and removed IdenticalIds and add the removed DNA node's own accession.
    new_ids = set(kept_ids) | set(removed_ids)
    new_ids.add(dna_node_to_remove["d.accession_id"])
    new_ids_list = list(new_ids)

    query_write_identical_ids = f"""
        MATCH (d:DNA) WHERE d.accession_id = "{dna_node_kept["d.accession_id"]}"
        SET d.IdenticalIds = {json.dumps(new_ids_list)}
    """
    eedb.db.execute_write(query_write_identical_ids)

    # Now remove the duplicate DNA node and reassign its relationships.
    remove_dna_node(dna_node_to_remove, dna_node_kept, eedb, logger)


def process_protein(current_protein_id):
    """
    Process a single protein by grouping its associated DNA nodes by the sliced sequence
    (from r.start to r.end) to quickly identify duplicates.
    """
    query_dna_nodes = f"""
        MATCH (reg:Region)-[rel_reg:HAS_REGION]-(d:DNA)-[r:ENCODES]->(p:Protein)
        WHERE p.accession_id = "{current_protein_id}" and reg.sequence_id = "{current_protein_id}"
        RETURN d.accession_id, d.sequence, rel_reg.start, rel_reg.end
    """
    dna_nodes = eedb.db.execute_read(query_dna_nodes)
    if not dna_nodes:
        return

    # Group DNA nodes by the sub-sequence slice.
    groups = {}
    for dna in dna_nodes:
        # Calculate the sub-sequence using r.start and r.end.
        sub_seq = dna["d.sequence"][dna["rel_reg.start"] : dna["rel_reg.end"]]
        groups.setdefault(sub_seq, []).append(dna)

    # For groups with duplicates, keep the first and remove the rest.
    for group in groups.values():
        if len(group) > 1:
            kept = group[0]
            for duplicate in group[1:]:
                # Skip if the accession IDs are the same
                if kept["d.accession_id"] == duplicate["d.accession_id"]:
                    continue

                LOGGER.info(
                    f"Identical DNA sequences: {kept['d.accession_id']} and {duplicate['d.accession_id']}"
                )
                update_and_remove_identical_dna(kept, duplicate, eedb, LOGGER)


def get_protein_ids_in_batches(batch_size=1000):
    """
    Generator that yields protein IDs in batches to avoid loading all protein IDs into memory.
    Uses SKIP and LIMIT in the Cypher query.
    """
    skip = 0
    while True:
        query = f"""
            MATCH (p:Protein)
            RETURN p.accession_id
            SKIP {skip} LIMIT {batch_size}
        """
        results = eedb.db.execute_read(query)
        if not results:
            break
        for record in results:
            yield record["p.accession_id"]
        skip += batch_size


if __name__ == "__main__":
    # Count of processed proteins.
    processed = 0
    # Process proteins in batches.
    for protein_id in get_protein_ids_in_batches(batch_size=200):
        process_protein(protein_id)
        processed += 1
        if processed % 100 == 0:
            LOGGER.info(f"Processed {processed} proteins.")

    # nohup python scr/code/011_TEM_Data_Cleaning_Smart_Combine_DNA.py > 011_TEM_Data_Cleaning_Smart_Combine_DNA.log 2>&1 &
