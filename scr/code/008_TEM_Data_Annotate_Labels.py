# ------------------------------------- IMPORTING PACKAGES -------------------------------------
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
    # read in the two dataframes form the blast searches
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
    # add to all of them the label "Source" with the value "BLAST_Protein"
    unique_ids_blast_protein = df_blast_protein["Subject ID"].unique()
    for id in unique_ids_blast_protein:
        # print(f"Annotating protein {id} with Source BLAST_Protein")
        label_node("Protein", id, "BLAST_Protein", "Source")
    LOGGER.info(
        f"Annotated {len(unique_ids_blast_protein)} proteins with Source BLAST_Protein"
    )

    # identify all unique ids in the df_tem_lactamase dataframe
    # add to all of them the label "Source" with the value "ELBD"
    unique_ids_elbd = df_tem_lactamase["protein_id_database"].unique()
    for id in unique_ids_elbd:
        # print(f"Annotating protein {id} with Source ELBD")
        label_node("Protein", id, "ELBD", "Source")
    LOGGER.info(f"Annotated {len(unique_ids_elbd)} proteins with Source ELBD")
    # identify all unique ids in the df_blast_dna dataframe
    # add to all of them the label "Source" with the value "BLAST_DNA"
    unique_ids_blast_dna = df_blast_dna["Subject ID"].unique()
    for id in unique_ids_blast_dna:
        # print(f"Annotating protein {id} with Source BLAST_DNA")
        label_node("DNA", id, "BLAST_DNA", "Source")
    LOGGER.info(f"Annotated {len(unique_ids_blast_dna)} proteins with Source BLAST_DNA")

    # identify all proteins that are not in the df_blast_protein dataframe or df_tem_lactamase dataframe
    # they where then added to the database via a DNA Connection
    # add to all of them the label "Source" with the value "DNA_Connection"

    # all proteins in the database
    query_all_proteins = "MATCH (p:Protein) RETURN p.accession_id"
    all_proteins = eedb.db.execute_read(query_all_proteins)
    all_proteins = [item["p.accession_id"] for item in all_proteins]

    # exclude the proteins that are in the df_blast_protein or df_tem_lactamase dataframe
    proteins_to_annotate = [
        id
        for id in all_proteins
        if id not in unique_ids_blast_protein and id not in unique_ids_elbd
    ]

    for id in proteins_to_annotate:
        # print(f"Annotating protein {id} with Source DNA_Connection")
        label_node("Protein", id, "DNA_Connection", "Source")

    LOGGER.info(
        f"Annotated {len(proteins_to_annotate)} proteins with Source DNA_Connection"
    )

    # identify all unique ids in the df_blast_dna dataframe
    # add to all of them the label "Source" with the value "BLAST_DNA"
    unique_ids_blast_dna = df_blast_dna["Subject ID"].unique()
    for id in unique_ids_blast_dna:
        # print(f"Annotating protein {id} with Source BLAST_DNA")
        label_node("DNA", id, "BLAST_DNA", "Source")
    LOGGER.info(f"Annotated {len(unique_ids_blast_dna)} proteins with Source BLAST_DNA")
