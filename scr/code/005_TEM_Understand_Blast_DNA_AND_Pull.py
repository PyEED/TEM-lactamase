# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os

import pandas as pd
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
    # df.to_csv(
    #     os.path.join(path_to_data_blast, "combined_data_blast_5000_tem_209.csv"),
    #     index=False,
    # )
    # read the csv file
    df = pd.read_csv(
        os.path.join(path_to_data_blast, "combined_data_blast_5000_dna_tem_209.csv")
    )
    print(df.head())
    print(len(df))
    # columns are: Query ID, Subject ID, E-value, Identity, Alignment Length, Mismatches, Gap Opens, Q. Start, Q. End, S. Start, S. End, E-value, Bit Score

    # how many unique subject id are there?
    unique_subject_ids = df["Subject ID"].unique()
    print(f"Number of unique subject ids: {len(unique_subject_ids)}")

    # check wether the ID is already in the database either in ids or in the IdenticalIds (if not Null)
    query_get_all_DNA_ids_and_identical_ids = """
        MATCH (d:DNA) RETURN d.accession_id, d.IdenticalIds
    """
    all_DNA_ids_and_identical_ids = eedb.db.execute_read(query_get_all_DNA_ids_and_identical_ids)
    print(f"Number of DNA entryies in the database: {len(all_DNA_ids_and_identical_ids)}")

    # create one long list of all the DNA ids including the IdenticalIds
    all_DNA_ids_plus_identical_ids = []
    for dna in all_DNA_ids_and_identical_ids:
        if dna["d.IdenticalIds"] is not None:
            all_DNA_ids_plus_identical_ids.extend(dna["d.IdenticalIds"])
        else:
            all_DNA_ids_plus_identical_ids.append(dna["d.accession_id"])
    print(f"Number of DNA ids plus identical ids: {len(all_DNA_ids_plus_identical_ids)}")

    # exclude the ones already in the database
    unique_subject_ids = [id for id in unique_subject_ids if id not in all_DNA_ids_plus_identical_ids]
    print(f"Number of unique subject ids after exclusion: {len(unique_subject_ids)}")

    for batch in range(0, len(unique_subject_ids), 500):
        trys = 0
        print(f"Batch {batch} of {len(unique_subject_ids)}")
        while trys < 20:
            try:
                batch_ids = unique_subject_ids[batch : batch + 500]

                eedb.fetch_from_primary_db(ids=batch_ids, db="ncbi_nucleotide")

                # if it run through the loop without an error, break the loop while loop
                break
            except Exception as e:
                LOGGER.error(f"Error processing batch {batch}: {str(e)}")
                trys += 1
                if trys == 20:
                    LOGGER.error("Maximum number of retries (20) reached")
                    break

    print("Done")


# nohup python scr/code/005_TEM_Understand_Blast_DNA_AND_Pull.py > 005_TEM_Understand_Blast_DNA_AND_Pull.log 2>&1 &
