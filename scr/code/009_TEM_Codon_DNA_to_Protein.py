# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os
import uuid

from dotenv import load_dotenv
from pyeed import Pyeed

# ------------------------------------- SETUP -------------------------------------


path_to_data_blast_dna = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-01-19_12-37-48"
path_to_data_blast_protein = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data/combined_data_blast_5000_tem_209"


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
        MATCH (r:Region {{annotation: "coding sequence"}})-[rel:HAS_REGION]-(d:DNA) WHERE d.accession_id = "{dna_id}" RETURN d.sequence, rel.start, rel.end
    """

    dna_sequence_query = eedb.db.execute_read(query_dna_sequence)
    if len(dna_sequence_query) == 0:
        LOGGER.info(f"No DNA sequence found for {dna_id}")
        return ""
    dna_sequence = dna_sequence_query[0]["d.sequence"]
    start = int(dna_sequence_query[0]["rel.start"])
    end = int(dna_sequence_query[0]["rel.end"])
    # Translate the DNA sequence using Biopython's Seq module.
    from Bio.Seq import Seq

    protein = str(Seq(dna_sequence[start:end]).translate(to_stop=True))
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

        # Check if this protein sequence already exists in the database
        query_check_existing_protein = f"""
            MATCH (p:Protein {{sequence: "{protein}"}})
            RETURN p.accession_id
        """
        existing_proteins = eedb.db.execute_read(query_check_existing_protein)

        if len(existing_proteins) > 0:
            # Protein sequence already exists, link DNA to existing protein
            existing_protein_id = existing_proteins[0]["p.accession_id"]
            LOGGER.info(f"Found existing protein with sequence: {existing_protein_id}")

            # Check if relationship already exists between DNA and protein
            query_check_relationship = f"""
                MATCH (d:DNA {{accession_id: "{dna_id}"}})-[r:ENCODES]->(p:Protein {{accession_id: "{existing_protein_id}"}})
                RETURN COUNT(r) AS relationship_count
            """
            relationship_exists = eedb.db.execute_read(query_check_relationship)

            if relationship_exists[0]["relationship_count"] == 0:
                # Create relationship between DNA and existing protein only if it doesn't exist
                query_dna_protein_relationship = f"""
                    MATCH (d:DNA) WHERE d.accession_id = "{dna_id}" 
                    MATCH (p:Protein) WHERE p.accession_id = "{existing_protein_id}" 
                    CREATE (d)-[:ENCODES {{start: 1, end: {len(protein)}}}]->(p)
                """
                eedb.db.execute_write(query_dna_protein_relationship)
                LOGGER.info(
                    f"Created relationship between DNA {dna_id} and existing protein {existing_protein_id}"
                )
            else:
                LOGGER.info(
                    f"Relationship already exists between DNA {dna_id} and protein {existing_protein_id}"
                )
        else:
            # Protein sequence doesn't exist, create new protein node
            id = uuid.uuid4()
            query_protein_node = f"""
                CREATE (p:Protein {{accession_id: "{id}", sequence: "{protein}", Source: "DE_NOVO_BASED_ON_DNA", seq_length: {len(protein)}}})
            """
            eedb.db.execute_write(query_protein_node)
            LOGGER.info(f"Created new protein node with accession id {id}")

            # Create relationship between DNA and new protein
            query_dna_protein_relationship = f"""
                MATCH (d:DNA) WHERE d.accession_id = "{dna_id}" 
                MATCH (p:Protein) WHERE p.accession_id = "{id}" 
                CREATE (d)-[:ENCODES {{start: 1, end: 1}}]->(p)
            """
            eedb.db.execute_write(query_dna_protein_relationship)
            LOGGER.info(
                f"Created relationship between DNA {dna_id} and new protein {id}"
            )


# nohup python scr/code/009_TEM_Codon_DNA_to_Protein.py > output_codon_dna_to_protein.log 2>&1 &
