# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os
import pandas as pd
from datetime import datetime

from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.tools.blast import Blast
from pyeed.analysis.standard_numbering import StandardNumberingTool
from pyeed.analysis.mutation_detection import MutationDetection


# ------------------------------------- SETUP -------------------------------------

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)


uri = "bolt://129.69.129.130:7687"
user = "neo4j"
eedb = Pyeed(uri, user=user, password='12345678')

eedb.calculate_sequence_embeddings(model_name="esmc_300m", batch_size=250)

"""
eedb.db.initialize_db_constraints(user, '12345678')

df = pd.read_csv(
    "/home/nab/Niklas/TEM-lactamase/data/002_combined_data/TEM_lactamase_with_dna_accession_id.csv", sep=","
)

df_protein_id_database = eedb.fetch_from_primary_db(
    df.loc[df['protein_id_database'].notna(), 'protein_id_database'].tolist(),
    db="ncbi_protein"
)

eedb.fetch_dna_entries_for_proteins()

# ------------------------------------- FUNCTIONS -------------------------------------

sn_protein = StandardNumberingTool(name="test_standard_numbering_protein")

sn_protein.apply_standard_numbering(
    base_sequence_id="AAP20891.1", db=eedb.db, list_of_seq_ids=df['protein_id_database'].tolist()
)

md = MutationDetection()

already_seen_pairs = set()

for protein_id_1 in df['protein_id_database'].dropna():
    for protein_id_2 in df['protein_id_database'].dropna():
        if protein_id_1 == protein_id_2:
            continue

        if (protein_id_1, protein_id_2) in already_seen_pairs:
            continue
        if (protein_id_2, protein_id_1) in already_seen_pairs:
            continue

        mutations_protein = md.get_mutations_between_sequences(
            protein_id_1, protein_id_2, eedb.db, sn_protein.name
        )

        already_seen_pairs.add((protein_id_1, protein_id_2))
        

"""