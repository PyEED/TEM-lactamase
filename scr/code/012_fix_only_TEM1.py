import logging
import os

from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.mutation_detection import MutationDetection
from pyeed.analysis.standard_numbering import StandardNumberingTool

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Base sequence (TEM-1a)
BASE_SEQUENCE_ID = "CAD09800.1"

# Name of the standard numbering tool that should already contain
# pairwise standard numbering between the base sequence and the flagged proteins.
# If the tool does not yet exist, it will be (re-)created automatically.
STANDARD_NUMBERING_TOOL_NAME = "standard_numbering_pairwise_flagged_proteins"

# Neo4j connection details are read from the environment (.env file)
load_dotenv()
PASSWORD = os.getenv("NEO4J_NIKLAS_TEM_THREE")
if PASSWORD is None:
    raise ValueError("NEO4J_NIKLAS_TEM_THREE is not set in the .env file.")

URI = "bolt://129.69.129.130:2137"
USER = "neo4j"

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def setup_logging() -> logging.Logger:
    """Configure root logger and return it."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    return logging.getLogger(__name__)


def get_flagged_proteins(db) -> list[str]:
    """Return a list of accession ids for proteins with ``length_flag = 1``."""
    query = (
        "MATCH (p:Protein)\n"
        "WHERE p.length_flag = 1\n"
        "RETURN p.accession_id AS protein_id"
    )
    results = db.execute_read(query)
    return [record["protein_id"] for record in results]


def ensure_standard_numbering(base_id: str, seq_ids: list[str], db) -> None:
    """Create or update pairwise standard numbering between *base_id* and *seq_ids*."""
    sn_tool = StandardNumberingTool(name=STANDARD_NUMBERING_TOOL_NAME)
    # The underlying implementation will silently ignore existing numbering, so we
    # can call this unconditionally.
    sn_tool.apply_standard_numbering_pairwise(
        base_sequence_id=base_id,
        db=db,
        list_of_seq_ids=seq_ids,
    )


def create_mutation_edges(base_id: str, target_ids: list[str], db, logger: logging.Logger) -> None:
    """Run mutation detection from *base_id* to every id in *target_ids*."""
    md_tool = MutationDetection()
    processed = 0
    for target_id in target_ids:
        if target_id == base_id:
            continue  # skip self comparison
        try:
            mutations = md_tool.get_mutations_between_sequences(
                sequence_id1=base_id,
                sequence_id2=target_id,
                db=db,
                standard_numbering_tool_name=STANDARD_NUMBERING_TOOL_NAME,
                save_to_db=True,
            )
            processed += 1
            logger.info(
                "Processed %s vs %s -> %d mutations",
                base_id,
                target_id,
                len(mutations["from_positions"]),
            )
        except Exception as exc:
            logger.error("Error processing %s vs %s: %s", base_id, target_id, exc)

    logger.info("Completed mutation detection for %d target proteins", processed)


# -----------------------------------------------------------------------------
# Main execution
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    logger = setup_logging()
    logger.info("Connecting to Neo4j at %s", URI)

    eedb = Pyeed(URI, user=USER, password=PASSWORD)

    try:
        # 1. obtain flagged proteins
        flagged_ids = get_flagged_proteins(eedb.db)
        logger.info("Retrieved %d length-flagged proteins", len(flagged_ids))

        # 2. ensure standard numbering exists
        ensure_standard_numbering(BASE_SEQUENCE_ID, flagged_ids, eedb.db)
        logger.info("Standard numbering ensured for base sequence against flagged proteins")

        # 3. run mutation detection
        create_mutation_edges(BASE_SEQUENCE_ID, flagged_ids, eedb.db, logger)
        logger.info("Mutation detection finished successfully")

    finally:
        # Always close the DB connection cleanly
        try:
            eedb.db.close()
        except Exception:
            pass 