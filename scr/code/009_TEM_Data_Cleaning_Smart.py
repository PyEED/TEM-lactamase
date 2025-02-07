# ------------------------------------- IMPORTING PACKAGES -------------------------------------
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
    Get all relationships for a node. So they they can be recreated for the other node that is kept.
    """
    query = f"""
    MATCH (n:{node_type} {{accession_id: "{node_id}"}})-[r]-() RETURN r
    """

    return eedb.db.execute_read(query)


if __name__ == "__main__":
    # here we are intrested in starting a data cleaning

    # first we want to find identical portein sequences

    # we get all of the portein_ids
    query_protein_ids = """
        MATCH (p:Protein) WHERE p.Source = "DE_NOVO_BASED_ON_DNA" RETURN p.accession_id
    """
    protein_ids = eedb.db.execute_read(query_protein_ids)
    protein_ids = [protein_id["p.accession_id"] for protein_id in protein_ids]
    print(f"Number of proteins: {len(protein_ids)}")

    et.drop_vector_index(index_name="vector_index_Protein_embedding", db=eedb.db)

    et.create_embedding_vector_index_neo4j(
        index_name="vector_index_Protein_embedding",
        db=eedb.db,
        similarity_function="cosine",
        m=512,
        ef_construction=3200,
        dimensions=960,
    )

    for index in range(0, len(protein_ids)):
        # print(f"Index {index} protein id: {protein_ids[index]}")

        results = et.find_nearest_neighbors_based_on_vector_index(
            index_name="vector_index_Protein_embedding",
            query_protein_id=protein_ids[index],
            number_of_neighbors=10,
            db=eedb.db,
        )
        if protein_ids[index] == "81bedbbf-597c-400a-a004-564f01ed0c9d":
            LOGGER.info(f"Index {index} protein id: {protein_ids[index]}")
            LOGGER.info(
                f"Found {len(results)} results for {protein_ids[index]} first match with score {results[0][1]}"
            )
            LOGGER.info(results)
        for result in results:
            if result[1] == 1.0:
                LOGGER.info(f"Found one that is identical to {protein_ids[index]}")

                # we now need to cleverly clean the data, since we do not want to lose any relationships

                # if one of the proteins is a DE_NOVO_BASED_ON_DNA then we want to remove that one and keep the other one
                # if both are DE_NOVO_BASED_ON_DNA then we just keep one of them

                current_protein_that_was_Searche_for_id = protein_ids[index]
                identical_protein_id = result[0]

                # we now need to check if the identical protein is a DE_NOVO_BASED_ON_DNA
                query_identical_protein = f"""
                    MATCH (p:Protein) WHERE p.accession_id = "{identical_protein_id}" RETURN p.Source
                """
                identical_protein_source = eedb.db.execute_read(query_identical_protein)

                if identical_protein_source[0]["p.Source"] == "DE_NOVO_BASED_ON_DNA":
                    identical_protein_is_de_novo = True
                else:
                    identical_protein_is_de_novo = False

                # we now need to check if the current protein is a DE_NOVO_BASED_ON_DNA
                query_current_protein = f"""
                    MATCH (p:Protein) WHERE p.accession_id = "{current_protein_that_was_Searche_for_id}" RETURN p.Source
                """
                current_protein_source = eedb.db.execute_read(query_current_protein)

                if current_protein_source[0]["p.Source"] == "DE_NOVO_BASED_ON_DNA":
                    current_protein_is_de_novo = True
                else:
                    current_protein_is_de_novo = False

                if identical_protein_is_de_novo and current_protein_is_de_novo:
                    LOGGER.info(
                        f"Both {current_protein_that_was_Searche_for_id} and {identical_protein_id} are DE_NOVO_BASED_ON_DNA"
                    )
                    # we now need to get all of the relationships for the identical protein that is then removed
                    identical_protein_relationships = get_all_realtiontship_for_node(
                        identical_protein_id, "Protein"
                    )
                    print(identical_protein_relationships)

                    # we now need to remove the identical protein
                    #
                    # query_remove_identical_protein = f"""
                    #     MATCH (p:Protein) WHERE p.accession_id = "{identical_protein_id}" DELETE p
                    # """
                    # eedb.db.execute_write(query_remove_identical_protein)

                    # we now need to recreate the relationships for the current protein
                    for relationship in identical_protein_relationships:
                        query_recreate_relationship = f"""
                            MATCH (p:Protein) WHERE p.accession_id = "{current_protein_that_was_Searche_for_id}"
                            CREATE (p)-[:{relationship["r"].type}]->({relationship["r"].end_node.id})
                        """
                        # eedb.db.execute_write(query_recreate_relationship)

                elif identical_protein_is_de_novo and not current_protein_is_de_novo:
                    LOGGER.info(
                        f"{current_protein_that_was_Searche_for_id} is not a DE_NOVO_BASED_ON_DNA and {identical_protein_id} is a DE_NOVO_BASED_ON_DNA"
                    )
                else:
                    LOGGER.info(
                        f"{current_protein_that_was_Searche_for_id} is a DE_NOVO_BASED_ON_DNA and {identical_protein_id} is not a DE_NOVO_BASED_ON_DNA"
                    )
