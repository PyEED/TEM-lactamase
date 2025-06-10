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
from pyeed.embeddings.processor import get_processor
import torch


# ------------------------------------- SETUP -------------------------------------

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)


uri = "bolt://129.69.129.130:7687"
user = "neo4j"
eedb = Pyeed(uri, user=user, password='12345678')

# eedb.db.initialize_db_constraints(user, '12345678')

#protein_name,phenotype,protein_id,protein_id_database,dna_accession_id
df = pd.read_csv(
    "/home/nab/Niklas/TEM-lactamase/data/002_combined_data/TEM_lactamase_with_dna_accession_id.csv", sep=","
)

# df_protein_id_database = eedb.fetch_from_primary_db(
#     df.loc[df['protein_id_database'].notna(), 'protein_id_database'].tolist(),
#     db="ncbi_protein"
# )

# i want a csv with the following columns:
# accession_id, sequence, seq_length, embedding, name

# Define the list of models
model_name_list = ['prot_t5_xl_bfd', "esmc_600m", "esmc_300m", "facebook/esm2_t33_650M_UR50D", 
                   'prot_t5_xl_uniref50', 'facebook/esm2_t36_3B_UR50D', 'facebook/esm2_t12_35M_UR50D', 
                   'facebook/esm2_t6_8M_UR50D', 'facebook/esm2_t30_150M_UR50D']

# Get the processor for embeddings
processor = get_processor()

# write the phenotype in the protein node of the db
query = """
MATCH (p:Protein)
WHERE p.accession_id = $accession_id
SET p.phenotype = $phenotype
"""
# Filter out any null/NaN values and only process rows where both fields exist
valid_data = df[df['protein_id_database'].notna() & df['phenotype'].notna()]
for _, row in valid_data.iterrows():
    accession_id = row['protein_id_database']
    phenotype = row['phenotype']
    eedb.db.execute_write(query, parameters={"accession_id": accession_id, "phenotype": phenotype})

# First get all proteins from the database
query = """
MATCH (p:Protein)
RETURN p.accession_id, p.sequence, p.seq_length, p.name, p.phenotype
"""

# run the query
result = eedb.db.execute_read(query)

# process results
result_df = pd.DataFrame(result, columns=["p.accession_id", "p.sequence", "p.seq_length", "p.name", "p.phenotype"])

# Initialize dictionary to store embeddings for each model
embeddings_dict = {}

# Process each model
for model_name in model_name_list:
    LOGGER.info(f"Processing model: {model_name}")
    
    # Get the model
    model = processor.get_or_create_model(model_name=model_name, device=torch.device('cuda:1'))
    
    # Initialize list to store embeddings for this model
    model_embeddings = []
    
    # Calculate embeddings for each sequence
    for seq in result_df['p.sequence']:
        last_layer_emb = model.get_single_embedding_last_hidden_state(sequence=seq)
        # Mean pool over the sequence length dimension
        mean_pooled = last_layer_emb.mean(axis=0)
        model_embeddings.append(mean_pooled)
    
    # Store embeddings in dictionary
    embeddings_dict[model_name] = model_embeddings
    
    # Clean up the model
    processor.remove_model(model_name)
    del model

# Add embeddings to the dataframe
for model_name, embeddings in embeddings_dict.items():
    result_df[f'embedding_{model_name}'] = embeddings

# save the result to a csv
result_df.to_csv("isaac_db_playground_read_out_embeddings.csv", index=False)

# Final cleanup
processor.cleanup()