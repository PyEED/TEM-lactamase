import logging
import multiprocessing
import os
import signal
import sys
import time
from functools import partial

# Limit joblib to a single worker inside each process to avoid runaway joblib/loky forks that
# otherwise exhaust CPU and RAM when libraries inside Pyeed use joblib implicitly.
os.environ["JOBLIB_MAX_WORKERS"] = "1"
try:
    from joblib import parallel_config
    parallel_config(n_jobs=1, prefer="threads")  # ensure any joblib call is effectively serial
except Exception:
    # If joblib is not available here, the environment variable is still enough
    pass

from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.mutation_detection import MutationDetection
from pyeed.analysis.standard_numbering import StandardNumberingTool

# --- Basic Setup ---
load_dotenv()
password = os.getenv("NEO4J_NIKLAS_TEM_THREE")
if password is None:
    raise ValueError("KEY is not set in the .env file.")

uri = "bolt://129.69.129.130:2137"
user = "neo4j"

# --- Helper Functions for Parallelization ---


def setup_logging(process_num):
    """Setup process-specific logging"""
    logging.basicConfig(
        level=logging.INFO,
        format=f"%(asctime)s - Process {process_num} - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(f"013_dna_mutations_process_{process_num}.log"),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger(__name__)


def process_part1_batch(process_num, batch, base_sequence_id, base_region_id, name_of_standard_numbering_tool, event):
    """Worker function for part 1 to process a batch of DNA sequences."""
    LOGGER = setup_logging(process_num)
    LOGGER.info(f"Starting Part 1 process {process_num} with {len(batch)} sequences")

    # Set up signal handler for this process
    def signal_handler(signum, frame):
        LOGGER.info(f"Process {process_num} received termination signal {signum}")
        event.set()
        sys.exit(0)

    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    
    if event.is_set():
        LOGGER.info(f"Process {process_num} termination event already set, exiting")
        return
        
    eedb = None
    try:
        eedb = Pyeed(uri, user=user, password=password)
        sn_dna = StandardNumberingTool(name=name_of_standard_numbering_tool)
        md = MutationDetection()

        region_ids = [item["elementId(r2)"] for item in batch]
        dna_accession_ids = [item["d2.accession_id"] for item in batch]
        
        LOGGER.info(f"Processing batch of {len(dna_accession_ids)} sequences.")

        sn_dna.apply_standard_numbering_pairwise(
            base_sequence_id=base_sequence_id,
            db=eedb.db,
            list_of_seq_ids=dna_accession_ids,
            node_type="DNA",
            region_ids_neo4j=region_ids + [base_region_id],
        )

        for i, accession_id in enumerate(dna_accession_ids):
            if event.is_set():
                LOGGER.info("Termination signal received, stopping.")
                break
            try:
                md.get_mutations_between_sequences(
                    sequence_id1=base_sequence_id,
                    sequence_id2=accession_id,
                    db=eedb.db,
                    node_type="DNA",
                    standard_numbering_tool_name=name_of_standard_numbering_tool,
                    region_ids_neo4j=[region_ids[i], base_region_id],
                )
            except Exception as e:
                LOGGER.error(f"Error processing sequence {accession_id}: {str(e)}")
                continue
        
        LOGGER.info(f"Process {process_num} completed Part 1 successfully")
        
    except Exception as e:
        LOGGER.error(f"Fatal error in Part 1 process {process_num}: {str(e)}")
        raise
    finally:
        if eedb:
            try:
                eedb.db.close()
            except:
                pass


def process_part2_batch(process_num, skip_idx, batch_size, name_of_standard_numbering_tool, key_base_sequence_id, key_base_sequence_Region_id, event):
    """Worker function for part 2 to process a batch of DNA pairs."""
    LOGGER = setup_logging(process_num)
    LOGGER.info(f"Starting Part 2 process {process_num} with skip index {skip_idx}")

    # Set up signal handler for this process
    def signal_handler(signum, frame):
        LOGGER.info(f"Process {process_num} received termination signal {signum}")
        event.set()
        sys.exit(0)

    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    
    if event.is_set():
        LOGGER.info(f"Process {process_num} termination event already set, exiting")
        return

    eedb = None
    try:
        eedb = Pyeed(uri, user=user, password=password)
        sn_dna = StandardNumberingTool(name=name_of_standard_numbering_tool)
        md = MutationDetection()

        query_protein_dna = """
        MATCH (r1:Region)-[:HAS_REGION]-(d1:DNA)-[:ENCODES]-(p1:Protein)-[:MUTATION]-(p2:Protein)-[:ENCODES]-(d2:DNA)-[:HAS_REGION]-(r2:Region)
        WHERE r1.annotation = $region_annotation
          AND r2.annotation = $region_annotation
          AND r1.sequence_id = p1.accession_id
          AND r2.sequence_id = p2.accession_id
          AND NOT d1.accession_id IN $list_of_seq_ids_to_skip
          AND NOT d2.accession_id IN $list_of_seq_ids_to_skip
          AND NOT (r1)-[:MUTATION]-(r2)
        RETURN DISTINCT elementId(r2), d2.accession_id, elementId(r1), d1.accession_id
        SKIP $skip
        LIMIT $limit
        """
        
        LOGGER.info(f"Processing batch with skip index {skip_idx}.")
        
        result_protein_dna = eedb.db.execute_read(
            query_protein_dna,
            parameters={
                "region_annotation": "coding sequence",
                "skip": skip_idx,
                "limit": batch_size,
                "list_of_seq_ids_to_skip": ["KT404265.1", "KT411823.1"],
            },
        )

        if not result_protein_dna:
            LOGGER.info(f"No more results at skip index {skip_idx}")
            return

        region_ids_r2 = [rec["elementId(r2)"] for rec in result_protein_dna]
        region_ids_r1 = [rec["elementId(r1)"] for rec in result_protein_dna]
        dna_accession_ids_target = [rec["d2.accession_id"] for rec in result_protein_dna]
        dna_accession_ids_source = [rec["d1.accession_id"] for rec in result_protein_dna]

        all_dna_ids = list(set(dna_accession_ids_target + dna_accession_ids_source))

        sn_dna.apply_standard_numbering_pairwise(
            base_sequence_id=key_base_sequence_id,
            db=eedb.db,
            list_of_seq_ids=all_dna_ids + [key_base_sequence_id],
            node_type="DNA",
            region_ids_neo4j=region_ids_r1 + region_ids_r2 + [key_base_sequence_Region_id],
        )

        for record in result_protein_dna:
            if event.is_set():
                LOGGER.info("Termination signal received, stopping.")
                break
            try:
                md.get_mutations_between_sequences(
                    sequence_id1=record["d1.accession_id"],
                    sequence_id2=record["d2.accession_id"],
                    db=eedb.db,
                    node_type="DNA",
                    standard_numbering_tool_name=name_of_standard_numbering_tool,
                    region_ids_neo4j=[record["elementId(r1)"], record["elementId(r2)"]],
                )
            except Exception as e:
                LOGGER.error(f"Error processing pair {record['d1.accession_id']} and {record['d2.accession_id']}: {str(e)}")
                continue
        
        LOGGER.info(f"Process {process_num} completed Part 2 successfully")
        
    except Exception as e:
        LOGGER.error(f"Fatal error in Part 2 process {process_num}: {str(e)}")
        raise
    finally:
        if eedb:
            try:
                eedb.db.close()
            except:
                pass


def main():
    """Main function to run the mutation detection pipeline."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    LOGGER = logging.getLogger(__name__)

    eedb = Pyeed(uri, user=user, password=password)
    
    # --- Find blaTEM1a database ID ---
    ids_tem = {}
    base_url_tem_family_card = "http://purl.obolibrary.org/obo/ARO_3000014"
    query_family = f"MATCH (o:OntologyObject {{name: '{base_url_tem_family_card}'}})-[*1..1]-(n) RETURN n"
    result_family = eedb.db.execute_read(query_family)

    for single_tem in result_family:
        if single_tem["n"]["name"] == "http://purl.obolibrary.org/obo/ARO_3000078":
            continue
        tem_name = single_tem["n"]["label"]
        tem_url = single_tem["n"]["name"]
        ids_tem[tem_name] = {"tem_name": tem_name, "tem_url": tem_url}
        query_protein = f"MATCH (o:OntologyObject {{name: '{tem_url}'}})-[*1..1]-(n:Protein) RETURN n"
        result_protein = eedb.db.execute_read(query_protein)
        if result_protein:
            protein_node = result_protein[0]["n"]
            ids_tem[tem_name]["identical_ids"] = protein_node.get("IdenticalIds", [])
            ids_tem[tem_name]["accession_id"] = protein_node.get("accession_id")

    blaTEM1a_id = "CAD09800.1"
    blaTEM1a_database_id = None
    for tem_data in ids_tem.values():
        if tem_data.get("accession_id") == blaTEM1a_id:
            blaTEM1a_database_id = tem_data["accession_id"]
            break
        if blaTEM1a_id in tem_data.get("identical_ids", []):
            blaTEM1a_database_id = tem_data["accession_id"]
            break
    
    if not blaTEM1a_database_id:
        raise ValueError("Could not find database ID for blaTEM1a.")
    LOGGER.info(f"Found blaTEM1a_database_id: {blaTEM1a_database_id}")

    # --- Configuration ---
    key_base_sequence_id = "AL513383.1"
    key_base_sequence_Region_id = '4:d0fe8add-d91b-4dda-9dd3-7b5b2e4fbfb9:23235041'
    name_of_standard_numbering_tool = "standard_numbering_pairwise_blaTEM1a_DNA_region_23235041_old_id"
    num_processes = 8
    batch_size = 1000  # Much larger batches = fewer total processes

    # --- Parallel Execution Setup ---
    termination_event = multiprocessing.Event()
    active_processes = []  # Track all active processes for proper cleanup
    
    def signal_handler(signum, frame):
        LOGGER.info(f"Main process received termination signal {signum}")
        termination_event.set()
        
        # Give processes time to clean up
        time.sleep(2)
        
        # Forcefully terminate any remaining processes
        for p in active_processes:
            if p.is_alive():
                p.terminate()
                p.join(timeout=1)
                if p.is_alive():
                    p.kill()
        
        sys.exit(0)
        
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

    # --- Part 1: Mutations vs Base Sequence ---
    LOGGER.info("--- Starting Part 1: Processing mutations against base sequence ---")
    query_vs_base = """
    MATCH (r1:Region)-[:HAS_REGION]-(d1:DNA)-[:ENCODES]-(p1:Protein)-[:MUTATION]-(p2:Protein)-[:ENCODES]-(d2:DNA)-[:HAS_REGION]-(r2:Region)
    WHERE p1.accession_id = $accession_id AND elementId(r1) = $region_id
      AND r2.sequence_id = p2.accession_id AND r1.annotation = 'coding sequence' AND r2.annotation = 'coding sequence'
      AND NOT (r1)-[:MUTATION]-(r2)
    RETURN DISTINCT elementId(r2), d2.accession_id
    """
    regions_to_process = eedb.db.execute_read(
        query_vs_base,
        parameters={"accession_id": blaTEM1a_database_id, "region_id": key_base_sequence_Region_id},
    )

    if regions_to_process:
        part1_batches = [regions_to_process[i:i+batch_size] for i in range(0, len(regions_to_process), batch_size)]
        
        # Process Part 1 batches in chunks of num_processes
        for chunk_start in range(0, len(part1_batches), num_processes):
            chunk_end = min(chunk_start + num_processes, len(part1_batches))
            current_chunk = part1_batches[chunk_start:chunk_end]
            
            # Create and start processes for this chunk
            chunk_processes = []
            for i, batch in enumerate(current_chunk):
                if batch:  # Only create process if batch is not empty
                    p = multiprocessing.Process(
                        target=process_part1_batch,
                        args=(
                            chunk_start + i,
                            batch,
                            key_base_sequence_id,
                            key_base_sequence_Region_id,
                            name_of_standard_numbering_tool,
                            termination_event
                        )
                    )
                    chunk_processes.append(p)
                    p.start()
                    active_processes.append(p) # Add to active processes
            
            try:
                # Wait for all processes in this chunk to complete
                for p in chunk_processes:
                    p.join()
                
                # Remove completed processes from active list
                for p in chunk_processes:
                    if p in active_processes:
                        active_processes.remove(p)
                        
            except Exception as e:
                LOGGER.error(f"Error during Part 1 processing chunk: {e}")
                termination_event.set()
                for p in chunk_processes:
                    p.join(timeout=5)
                    if p.is_alive():
                        p.terminate()
                        p.join(timeout=1)
                        if p.is_alive():
                            p.kill()
                    # Remove from active list
                    if p in active_processes:
                        active_processes.remove(p)
            
            LOGGER.info(f"Part 1 chunk {chunk_start//num_processes + 1} completed")
        
        LOGGER.info("Part 1 processing completed")
    else:
        LOGGER.info("No regions found for Part 1 processing.")

    # --- Part 2: All-vs-All Mutations ---
    LOGGER.info("--- Starting Part 2: Processing all-vs-all mutations ---")
    query_count = """
    MATCH (r1:Region)-[:HAS_REGION]-(d1:DNA)-[:ENCODES]-(p1:Protein)-[:MUTATION]-(p2:Protein)-[:ENCODES]-(d2:DNA)-[:HAS_REGION]-(r2:Region)
    WHERE r1.annotation = 'coding sequence' AND r2.annotation = 'coding sequence'
      AND NOT (r1)-[:MUTATION]-(r2)
    RETURN count(*) as total
    """
    # total_pairs_result = eedb.db.execute_read(query_count)
    # total_pairs = total_pairs_result[0]['total'] if total_pairs_result else 0
    total_pairs = 10_000_000  # Remove hardcoded value
    LOGGER.info(f"Total pairs: {total_pairs}")
    
    if total_pairs > 0:
        part2_work_items = [i for i in range(0, total_pairs, batch_size)]

        # Process Part 2 batches in chunks of num_processes
        for chunk_start in range(0, len(part2_work_items), num_processes):
            chunk_end = min(chunk_start + num_processes, len(part2_work_items))
            current_chunk = part2_work_items[chunk_start:chunk_end]
            
            # Create and start processes for this chunk
            chunk_processes = []
            for i, skip_idx in enumerate(current_chunk):
                p = multiprocessing.Process(
                    target=process_part2_batch,
                    args=(
                        chunk_start + i,
                        skip_idx,
                        batch_size,
                        name_of_standard_numbering_tool,
                        key_base_sequence_id,
                        key_base_sequence_Region_id,
                        termination_event
                    )
                )
                chunk_processes.append(p)
                p.start()
                active_processes.append(p) # Add to active processes

            try:
                # Wait for all processes in this chunk to complete
                for p in chunk_processes:
                    p.join()
                
                # Remove completed processes from active list
                for p in chunk_processes:
                    if p in active_processes:
                        active_processes.remove(p)
                        
            except Exception as e:
                LOGGER.error(f"Error during Part 2 processing chunk: {e}")
                termination_event.set()
                for p in chunk_processes:
                    p.join(timeout=5)
                    if p.is_alive():
                        p.terminate()
                        p.join(timeout=1)
                        if p.is_alive():
                            p.kill()
                    # Remove from active list
                    if p in active_processes:
                        active_processes.remove(p)
            
            LOGGER.info(f"Part 2 chunk {chunk_start//num_processes + 1} of {len(part2_work_items)//num_processes + 1} completed")

    else:
        LOGGER.info("No pairs found for Part 2 processing.")

    LOGGER.info("Done with all mutation detection tasks.")
    eedb.db.close()


if __name__ == "__main__":
    main()

# nohup python scr/code/013_TEM_Mutation_Detection_DNA.py > 013_TEM_Mutation_Detection_DNA.log 2>&1 &

# 2025-07-04 11:18:10.412 | DEBUG    | pyeed.analysis.mutation_detection:save_mutations_to_db:232 - Saved 6 mutations to database between AADNNK010000026.1 and KT398534.1
# MATCH (p1:Protein)-[:ENCODES]-(d1:DNA {accession_id: 'AADNNK010000026.1'})-[:HAS_REGION]-(r1:Region {annotation: 'coding sequence'})-[:MUTATION]-(r2:Region {annotation: 'coding sequence'})-[:HAS_REGION]-(d2:DNA {accession_id: 'KT398534.1'})-[:ENCODES]-(p2:Protein) RETURN p1, p2, d1, d2, r1, r2