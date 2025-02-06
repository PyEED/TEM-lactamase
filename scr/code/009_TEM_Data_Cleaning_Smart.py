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


if __name__ == "__main__":
    # here we are intrested in starting a data cleaning

    # first we want to find identical portein sequences

    # we get all of the portein_ids
    query_protein_ids = """
        MATCH (p:Protein) WHERE p.Source = "BLAST_Protein" RETURN p.accession_id
    """
    protein_ids = eedb.db.execute_read(query_protein_ids)
    protein_ids = [protein_id["p.accession_id"] for protein_id in protein_ids]
    print(f"Number of proteins: {len(protein_ids)}")

    et.drop_vector_index(index_name="vector_index_Protein_embedding", db=eedb.db)

    et.create_embedding_vector_index_neo4j(
        index_name="vector_index_Protein_embedding",
        db=eedb.db,
        similarity_function="cosine",
        m=48,
        ef_construction=300,
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
        if results[0][1] == 1.0:
            print(f"Found one that is identical to {protein_ids[index]}")
            print(results)
