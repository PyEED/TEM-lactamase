# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os

import pandas as pd
from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.ontology_loading import OntologyAdapter

# ------------------------------------- SETUP -------------------------------------


path_to_data_blast = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-01-19_12-37-48"
file_path = "/home/nab/Niklas/TEM-lactamase/CARD_Data_Ontologies/aro.owl"
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

# ------------------------------------- FUNCTIONS -------------------------------------


if __name__ == "__main__":
    # load in the the ontology of CARD, that means the structure of the ontology

    ontology_adaptor = OntologyAdapter()
    # ontology_adaptor.import_ontology_file_in_db(file_path, eedb.db)

    # now we want to pull all of the proteins that are in the CARD ontology, and link them to the ontology structure
    # we now open the tsv index file from CARD and link the proteins to the ontology, but first we have to pull them
    # ARO Accession	CVTERM ID	Model Sequence ID	Model ID	Model Name	ARO Name	Protein Accession	DNA Accession	AMR Gene Family	Drug Class	Resistance Mechanism	CARD Short Name

    # open the file and read in the proteins
    df = pd.read_csv(file_path_index, sep="\t")
    df = df.dropna(subset=["Protein Accession"])

    # now we want to fetch the proteins from the database
    eedb.fetch_from_primary_db(df["Protein Accession"].tolist(), db="ncbi_protein")

    # now we want to link the proteins to the ontology
    # we do this by matching the protein accession and the ARO Accession
    # the link realtionship is the following:     go_annotation = RelationshipTo("GOAnnotation", "ASSOCIATED_WITH")

    for index, row in df.iterrows():
        # the query is the following to match and to link the protein to the ontology
        # it is in cypther since we are using the neo4j database
        query = """
        MATCH (p:Protein {accession_id: $protein_accession})
        MATCH (a:OntologyObject {name: $aro_accession})
        MERGE (p)-[:ASSOCIATED_WITH]->(a)
        """

        # we now execute the query
        # example ARO:3002527 need to be name: http://purl.obolibrary.org/obo/ARO_0000001
        eedb.db.execute_write(
            query,
            parameters={
                "protein_accession": row["Protein Accession"],
                "aro_accession": f"http://purl.obolibrary.org/obo/{row['ARO Accession'].replace(':', '_')}",
            },
        )


# nohup python scr/code/001_TEM_Import_CARD.py > output_card_import.log 2>&1 &
