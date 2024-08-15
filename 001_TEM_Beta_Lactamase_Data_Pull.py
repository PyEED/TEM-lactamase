# this is the file that pulls all of our data

import os
import json
from pyeed.core import ProteinRecord

current_path = os.getcwd()

def read_in_proteins_from_text():
    # read in the text file
    file_path = '/home/nab/Niklas/TEM-lactamase/data/TEM_Ids/TEM_Ids.txt'
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

    print(protein_ids_fetched[0])

    return protein_ids_fetched

def remove_duplicates(proteins):
    
    # remove double entries
    data_protein_unique = []
    data_protein_unique_id = []
    for protein in proteins:
        if type(protein) != ProteinRecord:
            print(protein)
            continue
        if protein.id not in data_protein_unique_id:
            data_protein_unique.append(protein)
            data_protein_unique_id.append(protein.id)

    return data_protein_unique

def create_dic_for_protein_and_dna(proteins):

    # read in the file from json
    tem_dic = {}

    with open(os.path.join(current_path, "data", "TEM_Ids", "TEM_Ids.json"), "r") as f:
        tem_dic = json.load(f)

    data = {}

    # create a dic with one field protein and one field 'dna'
    for protein in proteins:
        try:
            data[tem_dic[protein.id.split('.')[0]]] = {
                'protein': protein,
                'dna': None
            }
        except KeyError:
            data[protein.id] = {
                'protein': protein,
                'dna': None
            }

    return data

def add_dna_data_in_dic(data):

    # now we pull the dna data
    for key in data.keys():
        protein = data[key]['protein']
        dna = protein.get_dna()
        data[key]['dna'] = dna

    return data

def save_output_to_json(data, folder_name):
    # save the blast search results
    output_folder_blast_search = os.path.join(current_path, "data", folder_name)
    os.makedirs(output_folder_blast_search, exist_ok=True)

    def dumper(obj):
        try:
            return obj.json()
        except:
            return obj.__dict__

    counter = 0
    for key, value in data.items():
        with open(os.path.join(output_folder_blast_search, key + ".json"), "w") as f:
            json.dump(value, f, default=dumper)
            counter += 1

def data_pull_all_run_local_blast_for_each_protein():

    proteins = read_in_proteins_from_text()
    proteins = remove_duplicates(proteins)
    data = create_dic_for_protein_and_dna(proteins)
    data = add_dna_data_in_dic(data)
    save_output_to_json(data, folder_name="TEM_Data_Pull_Only_IDS")
    
    results = []

    for protein in proteins:
        result = protein.ncbi_blast_local(db='/blast/blastdb/nr/nr', evalue = 0.005)
        results += result

        new_combined = result + proteins
        new_combined = remove_duplicates(new_combined)
        data = create_dic_for_protein_and_dna(new_combined)
        data = add_dna_data_in_dic(data)
        save_output_to_json(data, folder_name="TEM_Data_Pull_Latest_{}".format(protein.id))

if __name__ == "__main__":

    data_pull_all_run_local_blast_for_each_protein()