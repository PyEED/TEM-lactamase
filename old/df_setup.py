# read in the ids.json file form this directory
with open("/home/nab/Niklas/TEM-lactamase/data/TEM_Ids/TEM_Ids.json", "r") as f:
    dict_id_name = json.load(f)

print(dict_id_name)

# ----------------------------------------------------------------------------------------------------------------------------------------------

# Query to get all proteins from the tem id list with embeddings and get the label based on the annotations
# careful the dict has the key as AAP20891 and the database has a key as AAP20891.1
# we therefor load all accession_id from all proteins and then filter the proteins based on the accession_id
query = """
MATCH (p:Protein)
WHERE p.accession_id IS NOT NULL
RETURN p.accession_id AS protein_id
"""

results_all_ids = [record["protein_id"] for record in eedb.db.execute_read(query)]

ids = []

for i in results_all_ids:
    if i.split(".")[0] in dict_id_name.keys():
        ids.append(i)

print(f"Number of proteins in the database: {len(results_all_ids)}")
print(f"Number of proteins in the database with IDS matching the TEM ids: {len(ids)}")
print(f"Sample of ids: {ids[:5]}")


# ----------------------------------------------------------------------------------------------------------------------------------------------

file_path_label = '/home/nab/Niklas/TEM-lactamase/data/TEM_Ids/TEM_bldb_csv.csv'
df = pd.read_csv(file_path_label, sep=';')
# label columns is Bush_Jacoby_Class
id_labels = df[['Protein name', 'Phenotype']]
# convert it in a dict
dict_id_label = dict(zip(id_labels['Protein name'], id_labels['Phenotype']))
print(dict_id_label)

# the keys right now are TEM-1, TEM-2, etc. and the keys need to changed to AAP20891, from dict_id_name
dict_id_label_new = {}
for key in dict_id_name.keys():
    value = dict_id_name[key]
    # key is id, value is name
    if value in dict_id_label.keys():
        dict_id_label_new[key] = dict_id_label[value]

print(f"Sample of labels: {dict_id_label_new.keys()}")

# now there is the same problem a sbefore with the keys in the dict and the keys in the database
# we need to chnage the keys in the dict to match the keys in the database
dict_id_label_new_ids = {}
for i in results_all_ids:
    if i.split(".")[0] in dict_id_label_new.keys():
        dict_id_label_new_ids[i] = dict_id_label_new[i.split(".")[0]]

print(f"Sample of labels: {dict_id_label_new_ids.keys()}")
print(f"Sample of labels: {dict_id_label_new_ids}")

# some of them are still nan or ? both shoudl be None
dict_id_label_new_ids_none = {}
for key, value in dict_id_label_new_ids.items():
    if value == 'nan' or value == '?':
        dict_id_label_new_ids_none[key] = None
    # check nan but not as string, it is a pandas nan
    elif not pd.notna(value):
        dict_id_label_new_ids_none[key] = None
    else:
        dict_id_label_new_ids_none[key] = value

print(f"Sample of labels: {dict_id_label_new_ids_none}")

# next ensure that all of the proteins have a lable if it is missing do none
# dict_id_label_new_ids_none and ids have to have the same the missing on should be None
ids_with_labels = []
labels = []
for i in ids:
    if i not in dict_id_label_new_ids_none.keys():
        dict_id_label_new_ids_none[i] = None

print(f"Final length ids: {len(ids)} and labels {len(dict_id_label_new_ids_none)}")

# ----------------------------------------------------------------------------------------------------------------------------------------------

# all of the proteins need to be written into a pandas dataframe and then to a csv file
# the csv file should have the following columns: protein_id, protein_id_database, protein_name, phenotype
# the protein_id_database is the id from the database and the protein_id is the id from the dict
# the protein_name is the name from the dict and the phenotype is the label from the dict
# the csv file should be saved in the data folder

dataframe = []

# add protein_name and phenotype to the dataframe
for key in dict_id_label.keys():
    dataframe.append({'protein_name': key, 'phenotype': dict_id_label[key]})

print(dataframe[:5])

# add protein_id to the dataframe
for key in dict_id_name.keys():
    for row in dataframe:
        if row['protein_name'] == dict_id_name[key]:
            row['protein_id'] = key

print(dataframe[:5])

# add protein_id_database to the dataframe
for row in dataframe:
    for id in results_all_ids:
        if 'protein_id' in row.keys():
            if row['protein_id'] == id.split(".")[0]:
                row['protein_id_database'] = id


        

print(len(dataframe))

# write the dataframe to a csv file
df = pd.DataFrame(dataframe)
df.to_csv('/home/nab/Niklas/TEM-lactamase/data/002_combined_data/TEM_lactamase.csv', index=True, sep=';')

# ----------------------------------------------------------------------------------------------------------------------------------------------


