# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
from dotenv import load_dotenv
from pyeed import Pyeed

# ------------------------------------- SETUP -------------------------------------


path_to_data_blast = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-01-19_12-37-48"


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

# ------------------------------------- FUNCTIONS -------------------------------------


if __name__ == "__main__":
    # find out how many proteins are there which are connected to a DNA sequence
    query_ids_of_dna_connected_to_proteins = """
        MATCH (d:DNA)-[:ENCODES]->(p:Protein) RETURN d.accession_id
    """
    query_ids_of_dna_connected_to_proteins = eedb.db.execute_read(
        query_ids_of_dna_connected_to_proteins
    )
    print(
        f"Number of DNA sequences connected to a protein: {len(query_ids_of_dna_connected_to_proteins)}"
    )

    # find total number of proteins
    query_total_number_of_proteins = """
        MATCH (p:Protein) RETURN count(p)
    """
    total_number_of_proteins = eedb.db.execute_read(query_total_number_of_proteins)
    print(f"Total number of proteins: {total_number_of_proteins[0]['count(p)']}")

    # find total number of DNA sequences
    query_total_number_of_dna_sequences = """
        MATCH (d:DNA) RETURN count(d)
    """
    total_number_of_dna_sequences = eedb.db.execute_read(
        query_total_number_of_dna_sequences
    )
    print(
        f"Total number of DNA sequences: {total_number_of_dna_sequences[0]['count(d)']}"
    )

    # find the number of proteins which are standalone not conncted to DNA, mean no DNA connects to them
    query_number_of_proteins_without_dna_connecting_to_them = """
        MATCH (p:Protein) WHERE NOT (:DNA)-[:ENCODES]->(p) RETURN count(p)
    """
    number_of_proteins_without_dna_connecting_to_them = eedb.db.execute_read(
        query_number_of_proteins_without_dna_connecting_to_them
    )
    print(
        f"Number of proteins without DNA: {number_of_proteins_without_dna_connecting_to_them[0]['count(p)']}"
    )

    # find the number of DNA sequences which are not connected to a protein
    query_number_of_dna_sequences_without_protein_connecting_to_them = """
        MATCH (d:DNA) WHERE NOT (d)-[:ENCODES]->(:Protein) RETURN count(d)
    """
    number_of_dna_sequences_without_protein_connecting_to_them = eedb.db.execute_read(
        query_number_of_dna_sequences_without_protein_connecting_to_them
    )
    print(
        f"Number of DNA sequences without protein: {number_of_dna_sequences_without_protein_connecting_to_them[0]['count(d)']}"
    )

    # i want histogram frequency of the sequence length of the proteins, it is in the database as a property of the protein seq_length
    query_histogram_frequency_of_protein_sequence_length = """
        MATCH (p:Protein) RETURN p.seq_length
    """
    histogram_frequency_of_protein_sequence_length = eedb.db.execute_read(
        query_histogram_frequency_of_protein_sequence_length
    )
    # convert the list of dictionaries to a list of numpy arrays
    histogram_frequency_of_protein_sequence_length = [
        np.array(item["p.seq_length"])
        for item in histogram_frequency_of_protein_sequence_length
    ]

    # plot the histogram
    plt.hist(histogram_frequency_of_protein_sequence_length, bins=100)
    plt.title("Histogram of protein sequence length")
    plt.xlabel("Sequence length")
    plt.ylabel("Frequency")
    plt.savefig("histogram_frequency_of_protein_sequence_length.png")
    plt.close()

    # make second historgram just up the the seq length of 350
    plt.hist(histogram_frequency_of_protein_sequence_length, bins=100, range=(100, 450))
    plt.title("Histogram of protein sequence length up to 350")
    plt.xlabel("Sequence length")
    plt.ylabel("Frequency")
    plt.savefig("histogram_frequency_of_protein_sequence_length_up_to_350.png")
    plt.close()

    # a histogram which shows the number of connections a single DNA sequence has to other Proteins
    # to this for each DNA sequence in the database, count the number of proteins it connects to
    query_histogram_frequency_of_protein_connections_to_dna_sequence = """
        MATCH (d:DNA)
        OPTIONAL MATCH (d)-[:ENCODES]->(p:Protein)
        WITH d, COUNT(p) as num_connections
        RETURN num_connections, COUNT(d) as num_dnas
        ORDER BY num_connections
    """
    histogram_frequency_of_protein_connections_to_dna_sequence = eedb.db.execute_read(
        query_histogram_frequency_of_protein_connections_to_dna_sequence
    )

    # return the ids for 10 DNA sequences with the most connections to proteins
    query_ids_of_dna_sequences_with_most_connections_to_proteins = """
        MATCH (d:DNA)
        OPTIONAL MATCH (d)-[:ENCODES]->(p:Protein)
        WITH d, COUNT(p) as num_connections
        RETURN d.accession_id, num_connections
        ORDER BY num_connections DESC
        LIMIT 10
    """
    ids_of_dna_sequences_with_most_connections_to_proteins = eedb.db.execute_read(
        query_ids_of_dna_sequences_with_most_connections_to_proteins
    )
    print(
        f"IDs of DNA sequences with the most connections to proteins: {ids_of_dna_sequences_with_most_connections_to_proteins}"
    )

    # convert the list of dictionaries to a list of numpy arrays
    histogram_frequency_of_protein_connections_to_dna_sequence = [
        np.array(item["num_connections"])
        for item in histogram_frequency_of_protein_connections_to_dna_sequence
    ]

    # plot the histogram
    plt.hist(histogram_frequency_of_protein_connections_to_dna_sequence, bins=100)
    plt.title("Histogram of DNA sequences connected to proteins")
    plt.xlabel("Number of connections")
    plt.ylabel("Frequency")
    plt.savefig("histogram_frequency_of_dna_sequences_connected_to_proteins.png")
    plt.close()

    # plot once more for the number of connections up to 10, plot the absolute number of DNA sequences
    plt.bar(
        histogram_frequency_of_protein_connections_to_dna_sequence,
        histogram_frequency_of_protein_connections_to_dna_sequence,
    )
    plt.title("Histogram of DNA sequences connected to proteins up to 10 connections")
    plt.xlabel("Number of connections")
    plt.ylabel("Frequency")

# STAND 28.01.2025
# Number of DNA sequences connected to a protein: 24772
# Total number of proteins: 26118
# Total number of DNA sequences: 9981

# Number of DNA sequences connected to a protein: 205630
# Total number of proteins: 206976
# Total number of DNA sequences: 10910
