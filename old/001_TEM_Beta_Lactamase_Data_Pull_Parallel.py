# this file is used to pull the data from pyeed NCBI database
# the goal is to pull in parallel

# import the necessary libraries
import os
import sys
from pyeed.core import ProteinRecord
from concurrent.futures import ThreadPoolExecutor

sys.setrecursionlimit(250000)

# read in the text file
file_path = 'TEM_Ids.txt'
with open(file_path, 'r') as file:
    lines = file.readlines()

# create a dic with the name and id
# A	TEM-P118	t		AAN05029 	AY130285 	12354869 	view			2be	ESBL	A
# we want to keep the name and the id, these are the second and the fourth columns
# we need to split the line by the tab character
data = {}
count = 0
for line in lines:
    line = line.strip()
    columns = line.split('\t')
    name = columns[1]
    id = columns[4].strip()
    if id == 'EFJ83682':
        continue
    data[name] = id
    count += 1

protein_ids_fetched = ProteinRecord.get_ids(list(data.values()))

# here we want to add the region for the ncbi blast id
# this is a kind of work arround to get the id later for the name mappping
for i in range(0, len(protein_ids_fetched)):
    protein_ids_fetched[i].add_to_regions(
        name="blastquery:" + list(data.keys())[i],
        start=0,
        end=len(protein_ids_fetched[i].sequence)
    )


current_path = os.path.dirname(os.getcwd())
hits = 200

def fetch_protein(protein):
    try:
        return protein.ncbi_blast(n_hits=hits, e_value=0.1, db='nr')
    except Exception as e:
        print(f"Error fetching data for protein {protein.id}: {str(e)}")
        return None


with ThreadPoolExecutor(max_workers=10) as executor:
    results = list(executor.map(fetch_protein, protein_ids_fetched))


# Filter out None values which represent failed fetch attempts
results = [result for result in results if result is not None]


# add the protein_ids_fetched to results
for protein in protein_ids_fetched:
    results.append(protein)


# flatten the list if necessary and the item is from the List[ProteinRecord] type
data_proteins = []
for i in range(0, len(results)):
    if type(results[i]) == list:
        for protein in results[i]:
            if type(protein) == ProteinRecord:
                data_proteins.append(protein)
    else:
        if type(results[i]) == ProteinRecord:
            data_proteins.append(results[i])

# remove double entries
data_protein_unique = []
data_protein_unique_id = []
for protein in data_proteins:
    if type(protein) != ProteinRecord:
        print(protein)
        continue
    if protein.id not in data_protein_unique_id:
        data_protein_unique.append(protein)
        data_protein_unique_id.append(protein.id)



# save the blast search results
output_folder_blast_search = os.path.join(current_path, "TEM-lactamase", "data", "data_blast_search_nr_all_ids_hits_{}_db_{}".format(hits, len(data_protein_unique)))
os.makedirs(output_folder_blast_search, exist_ok=True)
counter = 0
for hit in data_proteins:
    with open(output_folder_blast_search + "/{}.json".format(counter), "w") as f:
        f.write(hit.json())
        counter += 1