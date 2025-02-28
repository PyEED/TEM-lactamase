# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os

import pandas as pd
from dotenv import load_dotenv
from pyeed import Pyeed

# ------------------------------------- SETUP -------------------------------------


path_to_data_blast_protein = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data/combined_data_blast_5000_tem_209"


load_dotenv()
password = os.getenv("NEO4J_NIKLAS_TEM_CLEAN")
if password is None:
    raise ValueError("KEY is not set in the .env file.")


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)


uri = "bolt://129.69.129.130:2123"
user = "neo4j"
eedb = Pyeed(uri, user=user, password=password)
eedb.db.initialize_db_constraints(user, password)

# ------------------------------------- FUNCTIONS -------------------------------------


def read_all_csv_files_and_concatenate_them(path_to_data_blast):
    # read all csv files in the directory and concatenate them
    all_files = [f for f in os.listdir(path_to_data_blast) if f.endswith(".csv")]
    df = pd.concat(
        [pd.read_csv(os.path.join(path_to_data_blast, f)) for f in all_files],
        ignore_index=True,
    )
    return df


if __name__ == "__main__":
    # df = read_all_csv_files_and_concatenate_them(path_to_data_blast)
    # save the dataframe to a csv file
    # df.to_csv(os.path.join(path_to_data_blast, 'combined_data_blast_5000_tem_209.csv'), index=False)
    # read the csv file
    df = pd.read_csv(
        os.path.join(path_to_data_blast_protein, "combined_data_blast_5000_tem_209.csv")
    )
    print(df.head())
    print(len(df))
    # columns are: Query ID, Subject ID, E-value, Identity, Alignment Length, Mismatches, Gap Opens, Q. Start, Q. End, S. Start, S. End, E-value, Bit Score

    # how many unique subject id are there?
    unique_subject_ids = df["Subject ID"].unique()
    print(f"Number of unique subject ids: {len(unique_subject_ids)}")

    for batch in range(0, len(unique_subject_ids), 2000):
        batch_ids = unique_subject_ids[batch : batch + 2000].tolist()
        eedb.fetch_from_primary_db(ids=batch_ids, db="ncbi_protein")
        eedb.fetch_dna_entries_for_proteins()


# nohup python scr/code/003_TEM_Understand_Blast_Protein_AND_Pull.py > output_protein_data_pull.log 2>&1 &
