# this is the big data pull for the beta lactamase data
# everything works now blast and backup so big stable data pull

# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os
from datetime import datetime

import numpy as np
import pandas as pd
from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.tools.blast import Blast

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


df = pd.read_csv(
    "/home/nab/Niklas/TEM-lactamase/data/002_combined_data/TEM_lactamase.csv", sep=";"
)
df["phenotype"] = df["phenotype"].apply(lambda x: np.nan if x == "?" else x)


path_to_data = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull"
path_to_data_blast = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data"
path_to_data_backup = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/backup_data"
path_to_data_log_file = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/log_file.txt"
path_to_db_blast = "/blast/db/custom/nr"

blast = Blast(
    mode="blastp",
    db_path=path_to_db_blast,
    db_name="nr",
    evalue=0.001,
    max_target_seqs=5000,
    num_threads=40,
)


def run_blast_and_fetch_data_proteins(
    id, path_to_data_blast, path_to_data_backup, path_to_data_log_file
):
    """
    This function runs the blast and fetches the data, in case of failure it will try again if that still does not work it will log and error but continue with the next id
    """
    eedb.fetch_from_primary_db(id, db="ncbi_protein")

    try:
        # query the seq of protein
        query = """
        MATCH (p:Protein)
        WHERE p.accession = $id
        RETURN p.sequence
        """

        # run the query
        sequence = eedb.db.execute_read(query, parameters={"id": id})

        # run the blast
        df_blast = blast.search(sequence)
        print(df_blast.head())
        df_blast.to_csv(f"{path_to_data_blast}/{id}.csv", index=False)
    except Exception as e:
        LOGGER.error(f"Error running blast for {id}: {e}")
        return None

    try:
        # try getting all of the ids in the neo4j
        # eedb.fetch_from_primary_db(df_blast['Subject ID'].tolist(), db='ncbi_protein')
        None

    except Exception as e:
        LOGGER.error(f"Error fetching data for {id}: {e}")
        return None


if __name__ == "__main__":
    # create folder of the timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    path_to_data_blast = f"{path_to_data_blast}/{timestamp}"
    path_to_data_backup = f"{path_to_data_backup}/{timestamp}"
    path_to_data_log_file = f"{path_to_data_log_file}/{timestamp}.txt"

    os.makedirs(path_to_data_blast, exist_ok=True)
    os.makedirs(path_to_data_backup, exist_ok=True)
    os.makedirs(path_to_data_log_file, exist_ok=True)

    for id in df["protein_id_database"].tolist():
        # checkout if the id is already been blasted
        path_previous_blast = f"/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data/2025-01-13_15-26-42/{id}.csv"
        if os.path.exists(path_previous_blast):
            print(f"Skipping {id} because it has already been blasted")
            continue

        # skip if is nan
        if pd.isna(id):
            continue
        print(str(id))
        print(type(id))
        run_blast_and_fetch_data_proteins(
            str(id), path_to_data_blast, path_to_data_backup, path_to_data_log_file
        )


# nohup python 002_TEM_Beta_Lactamase_Data_Pull.py > output.log 2>&1 &
