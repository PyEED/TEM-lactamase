# this is the big data pull for the beta lactamase data
# everything works now blast and backup so big stable data pull

# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os
import pandas as pd
from datetime import datetime

from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.tools.blast import Blast

# ------------------------------------- SETUP -------------------------------------

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)


uri = "bolt://129.69.129.130:7688"
user = "neo4j"
eedb = Pyeed(uri, user=user, password='12345678')
# eedb.db.wipe_database(date='2025-06-01')
# eedb.db.initialize_db_constraints(user, '12345678', models_path='/home/nab/Niklas/pyeed/src/pyeed/model.py')

path_to_data = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull"
path_to_data_blast_dna = (
    "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna"
)
path_to_data_backup = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/backup_data"
path_to_db_blast = "/databases/nt_new"
df = pd.read_csv(
    "/home/nab/Niklas/TEM-lactamase/data/002_combined_data/TEM_lactamase_with_dna_accession_id.csv", sep=","
)

blast = Blast(
    mode="blastn",
    db_path=path_to_db_blast,
    db_name="nt",
    evalue=0.001,
    max_target_seqs=5000,
    num_threads=30,
)

df_protein_id_database = eedb.fetch_from_primary_db(
    df.loc[df['protein_id_database'].notna(), 'protein_id_database'].tolist(),
    db="ncbi_protein"
)

eedb.fetch_dna_entries_for_proteins()

# ------------------------------------- FUNCTIONS -------------------------------------


def run_blast_and_fetch_data_dnas(id, path_to_data_blast_dna):
    """
    This function is supposed to identify all the dna sequences in the pull proteins in the database
    These DNA sequences are then blasted against the nt database. We want those blast results and save them to a csv file.
    """
    eedb.fetch_from_primary_db(id, db="ncbi_nucleotide")


    query_sequence = f"""
        MATCH (d:DNA) WHERE d.accession_id = "{id}"
        RETURN d.sequence
    """

    sequence = eedb.db.execute_read(query_sequence)[0]['d.sequence']
    LOGGER.info(f"Blasting {id}")
    df_blast = blast.search(sequence)

    LOGGER.info(f"Saving {id}")
    df_blast.to_csv(f"{path_to_data_blast_dna}/{id}.csv", index=False)


if __name__ == "__main__":

    # create folder of the timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    path_to_data_blast_dna = f"{path_to_data_blast_dna}/{timestamp}"
    path_to_data_backup = f"{path_to_data_backup}/{timestamp}"

    os.makedirs(path_to_data_blast_dna, exist_ok=True)
    os.makedirs(path_to_data_backup, exist_ok=True)

    query_ids_of_dna_connected_to_proteins = """
        MATCH (d:DNA)-[:ENCODES]->(p:Protein) RETURN d.accession_id
    """

    # /home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-05-16_14-03-46
    # check the ids from the folder
    ids_in_folder = os.listdir('/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-05-16_14-03-46')
    ids_in_folder = [id.split(".")[0] for id in ids_in_folder]
    ids_in_folder = [id + ".1" for id in ids_in_folder]
    print(ids_in_folder)


    query_ids_of_dna_connected_to_proteins = eedb.db.execute_read(
        query_ids_of_dna_connected_to_proteins
    )

    LOGGER.info(
        f"Number of DNA sequences to blast: {abs(len(query_ids_of_dna_connected_to_proteins)- len(ids_in_folder))}"
    )

    for id in query_ids_of_dna_connected_to_proteins:
        if id["d.accession_id"] in ids_in_folder:
            continue
        if id["d.accession_id"] == "LK391770.1":
            continue
        id_clean = id["d.accession_id"]
        print(f"Blasting {id_clean}")
        run_blast_and_fetch_data_dnas(id_clean, path_to_data_blast_dna)

    LOGGER.info("Done")
# nohup python scr/code/004_TEM_Beta_Lactamase_BLAST_DNA.py > 004_TEM_Beta_Lactamase_BLAST_DNA.log 2>&1 &
