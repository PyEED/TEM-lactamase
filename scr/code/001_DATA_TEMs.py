
import json
import logging
import pandas as pd
from pyeed import Pyeed

from pyeed.analysis.ontology_loading import OntologyAdapter
from pyeed.analysis.standard_numbering import StandardNumberingTool
from pyeed.analysis.sequence_alignment import PairwiseAligner

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
LOGGER = logging.getLogger(__name__)

# set up logger into a file
file_handler = logging.FileHandler("log.txt")

uri = "bolt://127.0.0.1:1123"
user = "neo4j"
password = "niklasonlytems"

# Create a Pyeed object, automatically connecting to the database
eedb = Pyeed(uri, user, password)
eedb.db.initialize_db_constraints(user=user, password=password)

# For testing purposes, we will wipe the database and remove all constraints
eedb.db.wipe_database(date='2022-12-32')
eedb.db.remove_db_constraints(user=user, password=password)

# DB connector is an attribute of the Pyeed object, type `DatabaseConnector`
LOGGER.info(f"Database stats: {eedb.db.stats()}")

# The first time the pyeed database is initialized, we need to create the constraints which are defined in the pyeed graph model
eedb.db.initialize_db_constraints(user=user, password=password)

# ok we are ready to go
LOGGER.info("Setup complete")

# read in the ids.json file form this directory
with open("/home/nab/Niklas/TEM-lactamase/data/TEM_Ids/TEM_Ids.json", "r") as f:
    dict_id_name = json.load(f)

# now fecth all of the proteins from the database
eedb.fetch_from_primary_db(dict_id_name, db='ncbi_protein')

# Apply the standard numbering
standard_numbering = StandardNumberingTool(name="test_standard_numbering")
standard_numbering.apply_standard_numbering(base_sequence_id='AAP20891.1', db=eedb.db)

# Align the sequences
aligner = PairwiseAligner()

# fetch all ids
query = """
        MATCH (p:Protein) 
        WHERE p.accession_id IS NOT NULL
        RETURN p.accession_id AS accession_id
        """
ids = [record['accession_id'] for record in eedb.db.execute_read(query)]
print(ids)

aligner.align_multipairwise(db=eedb.db, ids=ids)

# Fetch the DNA entries for the proteins
eedb.fetch_dna_entries_for_proteins()

# i want to know how many of the TEM proteins have a DNA sequence linked to them
# this can be found by checking if the DNA-[ENCODES]->Protein relationship exists
# then it should be compared to the TEM-Proteins from the dict and their IDs checked so that we can see if all of them have a DNA sequence

query = """
        MATCH (d:DNA)-[e:ENCODES]->(p:Protein)
        WHERE p.accession_id IS NOT NULL
        RETURN p.accession_id AS accession_id
        """

# we first need to compute the vector embedding for the proteins
eedb.calculate_sequence_embeddings(batch_size=32)

# loading in the owl ontology
file_path = "/home/nab/Niklas/TEM-lactamase/CARD_Data_Ontologies/aro.owl"
db = eedb.db
ontology_adapter = OntologyAdapter()

ontology_adapter.import_ontology_file_in_db(file_path, db)

# now we want to pull all of the proteins that are in the CARD ontology, and link them to the ontology structure
# we now open the tsv index file from CARD and link the proteins to the ontology, but first we have to pull them
# ARO Accession	CVTERM ID	Model Sequence ID	Model ID	Model Name	ARO Name	Protein Accession	DNA Accession	AMR Gene Family	Drug Class	Resistance Mechanism	CARD Short Name
file_path = "/home/nab/Niklas/TEM-lactamase/CARD_Data_Data/aro_index.tsv"

# open the file and read in the proteins
df = pd.read_csv(file_path, sep="\t")
df = df.dropna(subset=["Protein Accession"])

# now we want to fetch the proteins from the database
# eedb.fetch_from_primary_db(df["Protein Accession"].tolist(), db='ncbi_protein')

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
    # example ARO:3002527 need to be name: http://purl.obolibrary.org/obo/ARO_3002527
    eedb.db.execute_write(query, parameters={"protein_accession": row["Protein Accession"], "aro_accession": f"http://purl.obolibrary.org/obo/{row['ARO Accession'].replace(':', '_')}"})












