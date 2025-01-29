# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os

import numpy as np
from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.embedding_analysis import EmbeddingTool

# ------------------------------------- SETUP -------------------------------------


path_to_data_blast = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-01-19_12-37-48"


load_dotenv()
password = os.getenv("NEO4J_NIKLAS_TEM")
if password is None:
    raise ValueError("KEY is not set in the .env file.")


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)


uri = "bolt://127.0.0.1:8123"
user = "neo4j"
eedb = Pyeed(uri, user=user, password=password)
eedb.db.initialize_db_constraints(user, password)

et = EmbeddingTool()

# ------------------------------------- FUNCTIONS -------------------------------------


if __name__ == "__main__":
    # here we are intrested in starting a data cleaning

    # first we want to find identical portein sequences

    # we get all of the portein_ids
    query_protein_ids = """
        MATCH (p:Protein) RETURN p.accession_id
    """
    protein_ids = eedb.db.execute_read(query_protein_ids)
    print(f"Number of proteins: {len(protein_ids)}")

    index = 0

    print(f"Index {index} protein id: {protein_ids[index]}")

    # we check if the index protein has an embedding otherwise we keep searching for a index where there is an embedding
    query_embedding_exists = """
        MATCH (p:Protein {accession_id: $protein_id}) RETURN p.embedding IS NOT NULL
    """
    embedding_exists = eedb.db.execute_read(
        query_embedding_exists,
        parameters={"protein_id": protein_ids[index]["p.accession_id"]},
    )

    if embedding_exists[0]["p.embedding IS NOT NULL"]:
        print(f"Embedding exists for index {index}")
    else:
        print(f"Embedding does not exist for index {index}")

    # get the embedding of index as a numpy array
    query_embedding = """
        MATCH (p:Protein {accession_id: $protein_id}) RETURN p.embedding
    """
    embedding = eedb.db.execute_read(
        query_embedding,
        parameters={"protein_id": protein_ids[index]["p.accession_id"]},
    )

    embedding = np.array(embedding[0]["p.embedding"])
    print(f"Embedding: {embedding.shape}")

    results = et.find_closest_matches_simple(
        start_sequence_id=protein_ids[0]["p.accession_id"],
        db=eedb.db,
        metric="cosine",
        n=10,
    )

    print(results)
