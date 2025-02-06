# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os
import uuid

from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.embedding_analysis import EmbeddingTool
from pyeed.embedding import (
    load_model_and_tokenizer,
)

# ------------------------------------- SETUP -------------------------------------


path_to_data_blast_dna = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-01-19_12-37-48"
path_to_data_blast_protein = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data/combined_data_blast_5000_tem_209"


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

model, tokenizer, device = load_model_and_tokenizer("esmc_300m")


# ------------------------------------- FUNCTIONS -------------------------------------


def predict_protein_sequence(dna_id: str) -> str:
    """
    Predict the protein sequence encoded by the DNA id.
    This function looks for a FASTA file corresponding to the DNA accession id in the
    directory specified by 'path_to_data_blast_dna', reads the DNA sequence, and translates
    it into a protein sequence using the standard codon table.

    Parameters:
        dna_id (str): The accession id for the DNA sequence (expects a
                      file named '<dna_id>.fasta' within the specified directory).

    Returns:
        str: The translated protein sequence, or an empty string if the FASTA file is not found.
    """
    query_dna_sequence = f"""
        MATCH (d:DNA) WHERE d.accession_id = "{dna_id}" RETURN d.sequence
    """

    dna_sequence = eedb.db.execute_read(query_dna_sequence)
    dna_sequence = dna_sequence[0]["d.sequence"]
    # Translate the DNA sequence using Biopython's Seq module.
    from Bio.Seq import Seq

    protein = str(Seq(dna_sequence).translate(to_stop=True))
    return protein


# ------------------------------------- MAIN -------------------------------------
if __name__ == "__main__":
    # here we are intrested in starting a data cleaning
    # first we want to find identical portein sequences

    # we get all of the portein_ids
    query_dna_ids_without_protein = """
        MATCH (d:DNA) WHERE NOT (d)-[:ENCODES]->(:Protein) RETURN d.accession_id
    """
    dna_ids_without_protein = eedb.db.execute_read(query_dna_ids_without_protein)
    dna_ids_without_protein = [
        dna_id["d.accession_id"] for dna_id in dna_ids_without_protein
    ]
    print(f"Number of dna without protein: {len(dna_ids_without_protein)}")

    # create a for loop to predict the protein sequence for each dna id
    for dna_id in dna_ids_without_protein:
        LOGGER.info(f"Predicting protein sequence for {dna_id}")
        protein = predict_protein_sequence(dna_id)
        # convert the protein sequence to a string
        protein = str(protein)

        # we need to check if this protein is already in the database
        # beacuse the number of proteins is too large to check in a for loop we calulate the sequence embedding and check if it is already in the database
        # because we want to make our lives hard we first get the embedding, then save the protein to the database, then check it it already exisit, if it does then remove it
        # then either the newly added protein or the existing protein linked with DNA

        # now we want to generate the Protein Node with the newly generated protein sequence
        # the Source Attribute should be DE_NOVO_BASED_ON_DNA
        id = uuid.uuid4()
        query_protein_node = f"""
            CREATE (p:Protein {{accession_id: "{id}", sequence: "{protein}", Source: "DE_NOVO_BASED_ON_DNA", seq_length: {len(protein)}}})
        """
        eedb.db.execute_write(query_protein_node)

        eedb.calculate_sequence_embeddings(model_name="esmc_300m")

        # we need to check if the embedding is already in the database
        results = et.find_nearest_neighbors_based_on_vector_index(
            query_protein_id=id,
            index_name="vector_index_Protein_embedding",
            number_of_neighbors=2,
            db=eedb.db,
        )
        LOGGER.info(f"Results: {results}")

        already_in_database = False

        # if the score second attribute in the first return is 1.0 then we need to remove the protein
        if results[0][1] == 1.0:
            query_remove_protein = f"""
                MATCH (p:Protein) WHERE p.accession_id = "{id}" DELETE p
            """
            eedb.db.execute_write(query_remove_protein)
            LOGGER.info(f"Removed Protein Node with accession id {id}")
            already_in_database = True

        # now we want to create the relationship between the DNA and the Protein
        if already_in_database:
            query_dna_protein_relationship = f"""
                MATCH (d:DNA) WHERE d.accession_id = "{dna_id}" MATCH (p:Protein) WHERE p.accession_id = "{results[0][0]}" CREATE (d)-[:ENCODES]->(p)
            """
            LOGGER.info(
                f"Created relationship between DNA and Protein with accession id {results[0][0]} and DNA id {dna_id}"
            )
        else:
            query_dna_protein_relationship = f"""
                MATCH (d:DNA) WHERE d.accession_id = "{dna_id}" MATCH (p:Protein) WHERE p.accession_id = "{id}" CREATE (d)-[:ENCODES]->(p)
            """
            LOGGER.info(
                f"Created relationship between DNA and Protein with accession id {id} and DNA id {dna_id}"
            )

        eedb.db.execute_write(query_dna_protein_relationship)
