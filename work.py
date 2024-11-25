# this is the first general test for the fetching a protein and loading it into the database
import logging
from os import name
from turtle import st
from venv import logger

from pyeed.main import Pyeed
from pyeed.analysis.standard_numbering import StandardNumberingTool
from pyeed.analysis.sequence_alignment import PairwiseAligner

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
LOGGER = logging.getLogger(__name__)



def set_up_test():
    LOGGER.info("Setting up test")

    uri = "bolt://127.0.0.1:7687"
    user = "neo4j"
    password = "12345678"

    # Create a Pyeed object, automatically connecting to the database
    eedb = Pyeed(uri, user, password)

    # For testing purposes, we will wipe the database and remove all constraints
    eedb.db.wipe_database()
    eedb.db.remove_db_constraints(user=user, password=password)

    # DB connector is an attribute of the Pyeed object, type `DatabaseConnector`
    LOGGER.info(f"Database stats: {eedb.db.stats()}")

    # The first time the pyeed database is initialized, we need to create the constraints which are defined in the pyeed graph model
    eedb.db.initialize_db_constraints(user=user, password=password)

    return eedb

def run_work_flow(eedb):

    # ids 
    ids = [
        "AAL29438.1",
        'HBQ2613975.1',
        'EKW4005960.1',
        'EJG7116187.1',
        'AMM70781.1',
        'HCO3480053.1',
        'HAI5030310.1',
        'WP_000027057.1',
        'WP_215748091.1',
        'WP_261627585.1',
        'HBQ2613975.1',
        'EKW4005960.1',
        'EJG7116187.1',
        'AMM70781.1',
        'HCO3480053.1',
        'HAI5030310.1',
        'AII99784.1',
        'WP_000027057.1',
        'WP_215748091.1',
        'WP_261627585.1',
        'HDN1137928.1',
        'ARF47333.1',
        'WP_161654968.1',
        'EHC9934517.1',
        'ANG09566.1',
        'ANG13130.1',
        'ANG27879.1',
        'ARF46115.1',
        'EHI5535764.1',
        'WP_370620732.1',
        'ANG29652.1',
        'ANG25520.1',
        'ANG37455.1',
        'ANG13963.1',
        'ANG26903.1',
        'ANG18105.1',
        'ANG13537.1',
        'ANG30473.1',
        'ANG21110.1',
    ]
    


    # Fetch proteins from primary database
    eedb.fetch_from_primary_db(ids, db="ncbi_protein")

    # Check that the proteins were fetched
    LOGGER.info(f"Database stats: {eedb.db.stats()}")
    LOGGER.info("Test complete")

    ids_in_db = [i['p.accession_id'] for i in eedb.db.execute_read("MATCH (p:Protein) RETURN p.accession_id")]

    LOGGER.info(f"Proteins in database: {ids_in_db}")

    eedb.fetch_dna_entries_for_proteins()

    standard_numbering = StandardNumberingTool(name="test_standard_numbering")
    standard_numbering.apply_standard_numbering(base_sequence_id='AAL29438.1', db=eedb.db)

    # Align the sequences
    aligner = PairwiseAligner()
    aligner.align_multipairwise(db=eedb.db, ids=ids)


if __name__ == "__main__":
    eedb = set_up_test()

    # list all function in the eedb object
    LOGGER.info(dir(eedb))
    
    run_work_flow(eedb)

    # ids = ["NP_733140.2", "NP_733142.1", "VUD66739.1"]
    # eedb.fetch_from_primary_db(ids, db="NCBI")

    # ids = ["P04182","Q6QDP7"]

    # Fetch proteins from primary database
    # eedb.fetch_from_primary_db(ids)


    
    

    