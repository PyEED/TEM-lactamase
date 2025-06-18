import ast
import logging
import multiprocessing
import os
import signal
import sys
import time

import pandas as pd
from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.embedding_analysis import EmbeddingTool
from pyeed.analysis.mutation_detection import MutationDetection
from pyeed.analysis.standard_numbering import StandardNumberingTool

load_dotenv()
password = os.getenv("NEO4J_NIKLAS_TEM_THREE")
if password is None:
    raise ValueError("KEY is not set in the .env file.")

uri = "bolt://129.69.129.130:2137"
user = "neo4j"

blaTEM1a_database_id = "CAD09800.1"


def setup_logging(process_num):
    """Setup process-specific logging"""
    logging.basicConfig(
        level=logging.INFO,
        format=f"%(asctime)s - Process {process_num} - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(f"flagged_mutations_process_{process_num}.log"),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger(__name__)


def get_flagged_proteins():
    """Get all proteins with length_flag = 1"""
    eedb = Pyeed(uri, user=user, password=password)
    
    query = """
        MATCH (p:Protein)
        WHERE p.length_flag = 1
        RETURN p.accession_id as protein_id
    """
    
    results = eedb.db.execute_read(query)
    protein_ids = [record['protein_id'] for record in results]
    
    logging.info(f"Found {len(protein_ids)} flagged proteins")
    return protein_ids


def process_sequences(process_num, sequence_chunk, all_flagged_proteins, event):
    """Process a chunk of sequences in a single process"""

    # Set up signal handler for this process
    def signal_handler(signum, frame):
        LOGGER.info(f"Process {process_num} received termination signal {signum}")
        event.set()  # Set the event to notify other processes
        sys.exit(0)

    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

    LOGGER = setup_logging(process_num)
    LOGGER.info(f"Starting process {process_num} with {len(sequence_chunk)} sequences")

    eedb = Pyeed(uri, user=user, password=password)
    eedb.db.initialize_db_constraints(user, password)

    name_of_standard_numbering_tool = "standard_numbering_pairwise_flagged_proteins"
    et = EmbeddingTool()
    sn = StandardNumberingTool(name=name_of_standard_numbering_tool)
    md = MutationDetection()

    def process_mutation_pair(pair):
        try:
            mutations = md.get_mutations_between_sequences(
                pair[0],
                pair[1],
                eedb.db,
                name_of_standard_numbering_tool,
                save_to_db=True,
            )
            LOGGER.info(
                f"Processed pair {pair[0]} and {pair[1]}: {len(mutations['from_positions'])} mutations"
            )
            return mutations
        except Exception as e:
            LOGGER.error(f"Error processing pair {pair[0]} and {pair[1]}: {str(e)}")
            return None

    for count, id_base in enumerate(sequence_chunk):
        # Check if we should terminate
        if event.is_set():
            LOGGER.info(f"Process {process_num} received termination request")
            break

        LOGGER.info(
            f"Processing sequence {count + 1} of {len(sequence_chunk)} percentage {(count / len(sequence_chunk)) * 100:.2f}%"
        )

        # Find nearest neighbors for this protein
        try:
            results = et.find_nearest_neighbors_based_on_vector_index(
                eedb.db,
                id_base,
                index_name="vector_index_Protein_embedding",
                number_of_neighbors=len(all_flagged_proteins),
            )

            # Process mutations with other flagged proteins
            for result_id, result_distance in results:
                if result_id == id_base:
                    continue
                    
                # Only process if the result is also a flagged protein
                if result_id in all_flagged_proteins:
                    # To avoid processing the same pair twice, only process if id_base < result_id lexicographically
                    if id_base < result_id:
                        mutations = process_mutation_pair([id_base, result_id])
                        
        except Exception as e:
            LOGGER.error(f"Error processing sequence {id_base}: {str(e)}")
            continue

    LOGGER.info(f"Process {process_num} completed successfully")


def create_standard_numbering_tool(flagged_protein_ids):
    """Create standard numbering tool for flagged proteins"""
    eedb = Pyeed(uri, user=user, password=password)

    name_of_standard_numbering_tool = "standard_numbering_pairwise_flagged_proteins"
    sn = StandardNumberingTool(name=name_of_standard_numbering_tool)

    sn.apply_standard_numbering_pairwise(
        base_sequence_id=blaTEM1a_database_id, 
        db=eedb.db, 
        list_of_seq_ids=flagged_protein_ids
    )


if __name__ == "__main__":
    # Create an event to signal termination
    termination_event = multiprocessing.Event()

    def signal_handler(signum, frame):
        print(f"Main process received termination signal {signum}")
        termination_event.set()  # Set the event to notify all processes

        # Give processes time to clean up
        time.sleep(2)

        # Forcefully terminate any remaining processes
        for p in processes:
            if p.is_alive():
                p.terminate()
                p.join(timeout=1)
                if p.is_alive():
                    p.kill()  # Force kill if still alive

        sys.exit(0)

    # Set up signal handlers for the main process
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

    # Get flagged proteins from the database
    flagged_protein_ids = get_flagged_proteins()
    
    if not flagged_protein_ids:
        print("No flagged proteins found. Make sure to run the length flagging script first.")
        sys.exit(1)

    print(f"Found {len(flagged_protein_ids)} flagged proteins to process")

    # Create the standard numbering tool
    print("Creating standard numbering tool...")
    create_standard_numbering_tool(flagged_protein_ids)
    print("Standard numbering tool created")

    # Split sequences into 6 chunks
    chunk_size = len(flagged_protein_ids) // 6
    if chunk_size == 0:
        chunk_size = 1
        
    sequence_chunks = [
        flagged_protein_ids[i : i + chunk_size]
        for i in range(0, len(flagged_protein_ids), chunk_size)
    ]

    # Ensure we have at least one chunk
    if not sequence_chunks:
        sequence_chunks = [flagged_protein_ids]

    print(f"Split {len(flagged_protein_ids)} proteins into {len(sequence_chunks)} chunks")

    # Create and start processes
    processes = []
    for i, chunk in enumerate(sequence_chunks):
        if chunk:  # Only create process if chunk is not empty
            p = multiprocessing.Process(
                target=process_sequences, 
                args=(i, chunk, set(flagged_protein_ids), termination_event)
            )
            processes.append(p)
            p.start()

    try:
        # Wait for all processes to complete
        for p in processes:
            p.join()
    except Exception as e:
        print(f"Main process encountered an error: {e}")
        termination_event.set()
        # Wait for processes to terminate
        for p in processes:
            p.join(timeout=5)
            if p.is_alive():
                p.terminate()
                p.join(timeout=1)
                if p.is_alive():
                    p.kill()  # Force kill if still alive

    print("All processes completed successfully")

    # nohup python scr/code/012_Mutation_Detection_FLAGGED.py > 012_Mutation_Detection_FLAGGED.log 2>&1 &
