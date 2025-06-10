# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os
import uuid

from dotenv import load_dotenv
from pyeed import Pyeed

# ------------------------------------- SETUP -------------------------------------

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

# ------------------------------------- MAIN -------------------------------------

if __name__ == "__main__":
    # First, set all protein flags to 0
    query_reset_protein_flags = """
        MATCH (p:Protein)
        SET p.length_flag = 0
    """
    eedb.db.execute_write(query_reset_protein_flags)
    LOGGER.info("Reset all protein length flags to 0")

    # Set flag to 1 for proteins with correct length
    query_set_protein_flags = """
        MATCH (p:Protein) WHERE p.seq_length IN [285, 286, 287]
        SET p.length_flag = 1
        RETURN count(p) as updated_count
    """
    updated_proteins = eedb.db.execute_write(query_set_protein_flags)
    LOGGER.info(f"Set length_flag to 1 for {updated_proteins[0]['updated_count']} proteins")

    # First, set all DNA flags to 0
    query_reset_dna_flags = """
        MATCH (d:DNA)
        SET d.connected_flag = 0
    """
    eedb.db.execute_write(query_reset_dna_flags)
    LOGGER.info("Reset all DNA connected flags to 0")

    # Set flag to 1 for DNA that is connected to a protein with the correct length
    query_set_dna_flags = """
        MATCH (d:DNA)-[:ENCODES]->(p:Protein)
        WHERE p.seq_length IN [285, 286, 287]
        SET d.connected_flag = 1
        RETURN count(d) as updated_count
    """
    updated_dna = eedb.db.execute_write(query_set_dna_flags)
    LOGGER.info(f"Set connected_flag to 1 for {updated_dna[0]['updated_count']} DNA sequences")
