# this is the big data pull for the beta lactamase data
# everything works now blast and backup so big stable data pull

# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os
from datetime import datetime

from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.tools.blast import Blast

# ------------------------------------- SETUP -------------------------------------
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
blast = Blast()

path_to_data = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull"
path_to_data_blast_dna = (
    "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna"
)
path_to_data_backup = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/backup_data"
path_to_db_blast = "/blast/db/custom/nt_new/nt"

# ------------------------------------- FUNCTIONS -------------------------------------


def run_blast_and_fetch_data_dnas(id, path_to_data_blast_dna):
    """
    This function is supposed to identify all the dna sequences in the pull proteins in the database
    These DNA sequences are then blasted against the nt database. We want those blast results and save them to a csv file.
    """
    eedb.fetch_from_primary_db(id, db="ncbi_nucleotide")

    try:
        df_blast = blast.blastn(
            id,
            db=path_to_db_blast,
            dbConnector=eedb.db,
            evalue=0.001,
            outfmt=10,
            num_threads=40,
            max_target_seqs=5000,
        )
        print(df_blast.head())
        df_blast.to_csv(f"{path_to_data_blast_dna}/{id}.csv", index=False)
    except Exception as e:
        LOGGER.error(f"Error running blast for {id}: {e}")
        return None


if __name__ == "__main__":
    # eedb.fetch_dna_entries_for_proteins()

    # create folder of the timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    path_to_data_blast_dna = f"{path_to_data_blast_dna}/{timestamp}"
    path_to_data_backup = f"{path_to_data_backup}/{timestamp}"

    os.makedirs(path_to_data_blast_dna, exist_ok=True)
    os.makedirs(path_to_data_backup, exist_ok=True)

    query_ids_of_dna_connected_to_proteins = """
        MATCH (d:DNA)-[:ENCODES]->(p:Protein) RETURN d.accession_id
    """

    query_ids_of_dna_connected_to_proteins = eedb.db.execute_read(
        query_ids_of_dna_connected_to_proteins
    )

    print(
        f"Number of DNA sequences to blast: {len(query_ids_of_dna_connected_to_proteins)}"
    )

    for id in query_ids_of_dna_connected_to_proteins:
        id_clean = id["d.accession_id"]
        print(f"Blasting {id_clean}")
        run_blast_and_fetch_data_dnas(id_clean, path_to_data_blast_dna)


# nohup python 004_TEM_DNA_Blast.py > output.log 2>&1 &
