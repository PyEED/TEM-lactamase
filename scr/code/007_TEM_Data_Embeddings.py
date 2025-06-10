# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os

from dotenv import load_dotenv
from pyeed import Pyeed

# ------------------------------------- SETUP -------------------------------------


path_to_data_blast = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-01-19_12-37-48"


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

# ------------------------------------- FUNCTIONS -------------------------------------


if __name__ == "__main__":
    # calulcate the sequence embeddings on ems-c
    eedb.calculate_sequence_embeddings(model_name="esmc_300m", batch_size=250)


# nohup python scr/code/007_TEM_Data_Embeddings.py > output_embeddings.log 2>&1 &
