# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import concurrent.futures
import logging
import os

import pandas as pd
from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.embedding_analysis import EmbeddingTool

# ------------------------------------- SETUP -------------------------------------

path_to_data_blast_dna = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-01-19_12-37-48"
path_to_data_blast_protein = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data/combined_data_blast_5000_tem_209"
path_to_TEM_lactamase = (
    "/home/nab/Niklas/TEM-lactamase/data/002_combined_data/TEM_lactamase.csv"
)
file_path_index = "/home/nab/Niklas/TEM-lactamase/CARD_Data_Data/aro_index.tsv"


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


def parallel_label_nodes(node_type, ids, node_label, label_type, max_workers=10):
    """
    Process labeling tasks in parallel using ThreadPoolExecutor.

    Parameters:
      node_type (str): the type of the node, e.g. "Protein" or "DNA"
      ids (iterable): iterable of node IDs to annotate
      node_label (str): the label value to set
      label_type (str): the field/property to set on the node
      max_workers (int): number of concurrent workers.
    """
    futures_to_id = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        for _id in ids:
            futures_to_id[
                executor.submit(label_node, node_type, _id, node_label, label_type)
            ] = _id
        for future in concurrent.futures.as_completed(futures_to_id):
            _id = futures_to_id[future]
            try:
                future.result()
            except Exception as e:
                LOGGER.error(f"Error labeling {node_type} with id {_id}: {e}")


if __name__ == "__main__":
    # read in the two dataframes from the blast searches
    df_blast_dna = pd.read_csv(
        os.path.join(path_to_data_blast_dna, "combined_data_blast_5000_dna_tem_209.csv")
    )
    print(f"Head of df_blast_dna: {df_blast_dna.head()}")
    print(f"Shape of df_blast_dna: {df_blast_dna.shape}")
    df_blast_protein = pd.read_csv(
        os.path.join(path_to_data_blast_protein, "combined_data_blast_5000_tem_209.csv")
    )
    print(f"Head of df_blast_protein: {df_blast_protein.head()}")
    print(f"Shape of df_blast_protein: {df_blast_protein.shape}")
    # read in the TEM-lactamase dataframe from elbd
    df_tem_lactamase = pd.read_csv(path_to_TEM_lactamase, delimiter=";")
    print(f"Head of df_tem_lactamase: {df_tem_lactamase.head()}")
    print(f"Shape of df_tem_lactamase: {df_tem_lactamase.shape}")

    # identify all unique ids in the df_blast_protein dataframe
    unique_ids_blast_protein = df_blast_protein["Subject ID"].unique()
    # Annotate proteins with BLAST_Protein label in parallel
    parallel_label_nodes("Protein", unique_ids_blast_protein, "BLAST_Protein", "Source")
    LOGGER.info(
        f"Annotated {len(unique_ids_blast_protein)} proteins with Source BLAST_Protein"
    )

    # identify all the ids in the CARD index file
    df_card = pd.read_csv(file_path_index, sep="\t")
    df_card = df_card.dropna(subset=["Protein Accession"])
    unique_ids_card = df_card["Protein Accession"].unique()
    # Annotate proteins with CARD label in parallel
    parallel_label_nodes("Protein", unique_ids_card, "CARD", "Source")
    LOGGER.info(f"Annotated {len(unique_ids_card)} proteins with Source CARD")

    # identify all unique ids in the df_tem_lactamase dataframe
    unique_ids_elbd = df_tem_lactamase["protein_id_database"].unique()
    # Annotate proteins with BLDB label in parallel
    parallel_label_nodes("Protein", unique_ids_elbd, "BLDB", "Source")
    LOGGER.info(f"Annotated {len(unique_ids_elbd)} proteins with Source BLDB")

    # identify all unique ids in the df_blast_dna dataframe
    unique_ids_blast_dna = df_blast_dna["Subject ID"].unique()
    # Annotate DNA nodes with BLAST_DNA label in parallel
    parallel_label_nodes("DNA", unique_ids_blast_dna, "BLAST_DNA", "Source")
    LOGGER.info(
        f"Annotated {len(unique_ids_blast_dna)} DNA nodes with Source BLAST_DNA"
    )

    # identify all DNA which is not in the df_blast_dna dataframe, it then has to come from the Protein nucleotide ID
    # it should have the label DNA_FROM_PROTEIN
    query_all_dna = "MATCH (d:DNA) RETURN d.accession_id"
    all_dna = eedb.db.execute_read(query_all_dna)
    all_dna = [item["d.accession_id"] for item in all_dna]
    dna_to_annotate = [id for id in all_dna if id not in unique_ids_blast_dna]
    # Annotate these DNA nodes with BLAST_DNA label in parallel
    parallel_label_nodes("DNA", dna_to_annotate, "DNA_FROM_PROTEIN", "Source")
    LOGGER.info(
        f"Annotated {len(dna_to_annotate)} DNA nodes with Source DNA_FROM_PROTEIN"
    )

    # identify all proteins that are not in the df_blast_protein dataframe or df_tem_lactamase dataframe
    query_all_proteins = "MATCH (p:Protein) RETURN p.accession_id"
    all_proteins = eedb.db.execute_read(query_all_proteins)
    all_proteins = [item["p.accession_id"] for item in all_proteins]

    proteins_to_annotate = [
        id
        for id in all_proteins
        if id not in unique_ids_blast_protein
        and id not in unique_ids_elbd
        and id not in unique_ids_card
    ]
    # Annotate these proteins with PREDICTED_FROM_DNA label in parallel
    # because if not from ELDB not from CARD and not from BLAST it has to be predicted from DNA
    parallel_label_nodes(
        "Protein", proteins_to_annotate, "PREDICTED_FROM_DNA", "Source"
    )
    LOGGER.info(
        f"Annotated {len(proteins_to_annotate)} proteins with Source PREDICTED_FROM_DNA"
    )

    eedb.create_coding_sequences_regions()

# nohup python scr/code/008_TEM_Data_Annotate_Labels.py > output_annotate_labels.log 2>&1 &
