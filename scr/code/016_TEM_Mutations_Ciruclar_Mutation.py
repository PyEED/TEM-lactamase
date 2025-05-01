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

path_to_data_blast = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-01-19_12-37-48"

load_dotenv()
password = os.getenv("NEO4J_NIKLAS_TEM_NEW_START")
if password is None:
    raise ValueError("KEY is not set in the .env file.")

uri = "bolt://129.69.129.130:2127"
user = "neo4j"

blaTEM1a_database_id = "WP_000027057.1"


def setup_logging(process_num):
    """Setup process-specific logging"""
    logging.basicConfig(
        level=logging.INFO,
        format=f"%(asctime)s - Process {process_num} - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(f"process_{process_num}.log"),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger(__name__)


def process_sequences(process_num, sequence_chunk, event):
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

    name_of_standard_numbering_tool = (
        "standard_numbering_pairwise_circular_mutations_to_blaTEM1a"
    )
    et = EmbeddingTool()
    sn = StandardNumberingTool(name=name_of_standard_numbering_tool)
    md = MutationDetection()

    max_number_of_mutations = 10

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
            f"Processing sequence {count + 1} of {len(sequence_chunk)} percentage {(count / len(sequence_chunk)) * 100}"
        )
        found_more_than_max_number_of_mutations = False
        start_number = 1000
        count_break = 0

        while not found_more_than_max_number_of_mutations:
            results = et.find_nearest_neighbors_based_on_vector_index(
                eedb.db,
                id_base,
                index_name="vector_index_Protein_embedding",
                number_of_neighbors=start_number,
            )

            for result_id, result_distance in results:
                if result_id == id_base:
                    continue
                if result_id in sequence_chunk:
                    mutations = process_mutation_pair([id_base, result_id])
                    count_break += 1
                    if mutations is not None:
                        if len(mutations["from_positions"]) > max_number_of_mutations:
                            found_more_than_max_number_of_mutations = True
                            LOGGER.info(
                                f"Found more than {max_number_of_mutations} mutations for {id_base} and {result_id} after {count_break} neighbors"
                            )
                            break

            start_number += 5000

    LOGGER.info(f"Process {process_num} completed successfully")


def create_standard_numbering_tool(ids_in_circle):
    # connect to the database

    eedb = Pyeed(uri, user=user, password=password)

    name_of_standard_numbering_tool = (
        "standard_numbering_pairwise_circular_mutations_to_blaTEM1a"
    )
    sn = StandardNumberingTool(name=name_of_standard_numbering_tool)

    sn.apply_standard_numbering_pairwise(
        base_sequence_id=blaTEM1a_database_id, db=eedb.db, list_of_seq_ids=ids_in_circle
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

    # Read and prepare data
    path = "/home/nab/Niklas/TEM-lactamase/combined_protein_data.csv"
    df = pd.read_csv(path)

    names = [
        "TEM beta-lactamase",
        "GIL beta-lactamase",
        "AER beta-lactamase",
        "LAP beta-lactamase",
        "HERA beta-lactamase",
        "CKO beta-lactamase",
        "GES beta-lactamase",
    ]

    df_tem = df[df["family_name"].isin(names)]
    ids_in_circle_series = df_tem["ids_in_circle"].apply(
        lambda x: ast.literal_eval(x) if x != "[]" else []
    )
    ids_in_circle = []
    for circle_list in ids_in_circle_series:
        ids_in_circle.extend(circle_list)
    ids_in_circle = list(dict.fromkeys(ids_in_circle))

    # create the standard numbering tool
    create_standard_numbering_tool(ids_in_circle)

    # Split sequences into 6 chunks
    chunk_size = len(ids_in_circle) // 6
    sequence_chunks = [
        ids_in_circle[i : i + chunk_size]
        for i in range(0, len(ids_in_circle), chunk_size)
    ]

    # Ensure we have exactly 6 chunks (last chunk might be smaller)
    while len(sequence_chunks) < 6:
        sequence_chunks.append([])

    # Create and start processes
    processes = []
    for i in range(6):
        p = multiprocessing.Process(
            target=process_sequences, args=(i, sequence_chunks[i], termination_event)
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

    # nohup python scr/code/016_TEM_Mutations_Ciruclar_Mutation.py > 016_TEM_Mutations_Ciruclar_Mutation.log 2>&1 &
