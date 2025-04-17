import logging
import os

from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.embedding_analysis import EmbeddingTool

path_to_data_blast = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-01-19_12-37-48"


load_dotenv()
password = os.getenv("NEO4J_NIKLAS_TEM_CLEAN")
if password is None:
    raise ValueError("KEY is not set in the .env file.")


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)


uri = "bolt://129.69.129.130:2123"
user = "neo4j"
eedb = Pyeed(uri, user=user, password=password)
eedb.db.initialize_db_constraints(user, password)

# First, collect all TEM-lactamase proteins and their identical IDs
data = {}

# Base URLs for ontology queries
base_url_tem_family_card = "http://purl.obolibrary.org/obo/ARO_3000014"
base_url_class_A = "http://purl.obolibrary.org/obo/ARO_3000078"


# Helper function to handle family terms
def handle_family_term(single_family):
    query = f"""
    MATCH (o:OntologyObject {{name: '{single_family["n"]["name"]}'}})-[*1..1]-(n) RETURN n
    """
    result = eedb.db.execute_read(query)
    return result


# Now get all Class A beta-lactamase family proteins
query = f"""
MATCH (o:OntologyObject {{name: '{base_url_class_A}'}})-[*1..1]-(n) RETURN n
"""
result = eedb.db.execute_read(query)

for single_family in result:
    family_name = single_family["n"]["label"]
    family_result = handle_family_term(single_family)

    # Special handling for SHV-LEN family
    if single_family["n"]["name"] == "http://purl.obolibrary.org/obo/ARO_3000096":
        query = f"""
        MATCH (o:OntologyObject {{name: '{single_family["n"]["name"]}'}})-[*1..1]-(n) RETURN n
        """
        family_result = eedb.db.execute_read(query)

        family_result = handle_family_term(family_result[0])
        family_result = family_result + handle_family_term(family_result[1])

    for single_protein in family_result:
        protein_name = single_protein["n"]["label"]
        protein_url = single_protein["n"]["name"]

        query_protein_url = f"""
        MATCH (o:OntologyObject {{name: '{protein_url}'}})-[*1..1]-(n:Protein) RETURN n
        """
        result_protein_url = eedb.db.execute_read(query_protein_url)

        if len(result_protein_url) == 0:
            continue

        protein_id = result_protein_url[0]["n"]["accession_id"]

        data[protein_name] = {
            "protein_id": protein_id,
            "protein_name": protein_name,
            "protein_url": protein_url,
            "family_name": family_name,
        }

# Initialize embedding tool
et = EmbeddingTool()

# Create a working dataset with proportional representation of each family
data_working = {}
family_count = {}
total_proteins = 0

# Count proteins per family and total proteins
for protein_name, protein_data in data.items():
    family_name = protein_data["family_name"]
    if family_name not in family_count:
        family_count[family_name] = 0
    family_count[family_name] += 1
    total_proteins += 1

# Calculate the proportion for each family
# We want to keep a reasonable total number of proteins (around 100-150)
target_total = min(150, total_proteins)
family_max_proteins = {}

for family_name, count in family_count.items():
    # Ensure each family gets at least one protein
    # For larger families, allocate proportionally more proteins
    proportion = count / total_proteins
    max_proteins = min(3, max(1, round(proportion * target_total)))
    family_max_proteins[family_name] = max_proteins

# Add proteins to the working dictionary based on calculated proportions
for family_name in family_count:
    family_proteins_added = 0

    for protein_name, protein_data in data.items():
        if (
            protein_data["family_name"] == family_name
            and family_proteins_added < family_max_proteins[family_name]
        ):
            data_working[protein_name] = protein_data
            family_proteins_added += 1

print("Proteins per family:")
for family, max_count in family_max_proteins.items():
    print(f"{family}: {max_count}")
print(f"Total proteins after proportional selection: {len(data_working)}")

if __name__ == "__main__":
    query_number_of_total_ids = "MATCH (n:Protein) RETURN count(n)"
    result_number_of_total_ids = eedb.db.execute_read(query_number_of_total_ids)
    number_of_total_ids = result_number_of_total_ids[0]["count(n)"]

    # Initialize lists
    already_processed_ids = []

    minimal_cosine_similarity_in_circle = [1 for _ in range(len(data_working))]
    ids_in_circle = [[] for _ in range(len(data_working))]
    cosine_similarity_scores_vs_neighbors_data_lines = [
        [1] for _ in range(len(data_working))
    ]
    number_of_neighbours = [[1] for _ in range(len(data_working))]
    neighbour_added = [1 for _ in range(len(data_working))]

    # Create a stable mapping of protein names to indices
    protein_indices = {name: idx for idx, name in enumerate(data_working.keys())}
    first_round = True
    counter = 0
    while len(already_processed_ids) < int(number_of_total_ids * 0.8):
        counter += 1
        print(
            f"Processing batch {counter} total number of ids: {number_of_total_ids} currently processed: {len(already_processed_ids)}"
        )

        global_minal_this_round = min(minimal_cosine_similarity_in_circle) - 0.01

        # Process all proteins sequentially
        for protein_name, protein_data in data_working.items():
            index_in_list = protein_indices[protein_name]
            protein_id = protein_data["protein_id"]
            protein_url = protein_data["protein_url"]
            family_name = protein_data["family_name"]

            result = et.find_nearest_neighbors_based_on_vector_index(
                index_name="vector_index_Protein_embedding",
                db=eedb.db,
                query_protein_id=protein_id,
                number_of_neighbors=number_of_neighbours[index_in_list][-1]
                + neighbour_added[index_in_list],
            )

            neighbour_added[index_in_list] = 1

            timer = 0

            while global_minal_this_round <= result[-1][1] or first_round:
                # if neighbour_added[index_in_list] > 100 and counter < 3:
                #    break
                if first_round:
                    first_round = False

                print(
                    f"Neighbour search radius: {number_of_neighbours[index_in_list][-1] + neighbour_added[index_in_list]} for protein {protein_name} with id {protein_id} with minimal cosine similarity {global_minal_this_round} currennt result[-1]: {result[-1]}"
                )
                neighbour_no = 0
                found_minimum = False

                for i in result:
                    neighbour_no += 1
                    if i[0] not in already_processed_ids:
                        if i[1] >= global_minal_this_round or (
                            found_minimum and neighbour_no <= len(result)
                        ):
                            already_processed_ids.append(i[0])
                            ids_in_circle[index_in_list].append(i[0])
                            number_of_neighbours[index_in_list].append(neighbour_no)
                            cosine_similarity_scores_vs_neighbors_data_lines[
                                index_in_list
                            ].append(i[1])

                        if i[1] < global_minal_this_round:
                            found_minimum = True

                        if (
                            minimal_cosine_similarity_in_circle[index_in_list]
                            > cosine_similarity_scores_vs_neighbors_data_lines[
                                index_in_list
                            ][-1]
                        ):
                            minimal_cosine_similarity_in_circle[index_in_list] = (
                                cosine_similarity_scores_vs_neighbors_data_lines[
                                    index_in_list
                                ][-1]
                            )

                neighbour_added[index_in_list] += 1 + timer**2

                result = et.find_nearest_neighbors_based_on_vector_index(
                    index_name="vector_index_Protein_embedding",
                    db=eedb.db,
                    query_protein_id=protein_id,
                    number_of_neighbors=number_of_neighbours[index_in_list][-1]
                    + neighbour_added[index_in_list],
                )
                timer += 1

        # plot all the data
        print(f"Number of cirlces found: {len(ids_in_circle)}")
        print(f"Number of already processed ids: {len(already_processed_ids)}")
        print(f"Number of total ids: {number_of_total_ids}")
        print(f"Number of neighbours: {number_of_neighbours}")

        import matplotlib.pyplot as plt

        plt.figure(figsize=(10, 10))
        plt.plot(
            minimal_cosine_similarity_in_circle,
        )
        plt.savefig("minimal_cosine_similarity_in_circle.png")
        plt.close()
        # plot the neighbours vs the cosine similarity
        plt.figure(figsize=(15, 10))  # Increased width to accommodate legend

        for i, (neighbors, similarities) in enumerate(
            zip(number_of_neighbours, cosine_similarity_scores_vs_neighbors_data_lines)
        ):
            plt.plot(
                neighbors[1:],
                similarities[1:],
                label=f"{list(data_working.values())[i]['family_name']}",
            )

        plt.xlabel("Number of Neighbors")
        plt.ylabel("Cosine Similarity")
        plt.title("Neighbors vs Cosine Similarity")

        # Create legend outside the plot
        legend = plt.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

        # Adjust layout to make room for the legend
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # Left, bottom, right, top

        # Save the figure with the legend included
        plt.savefig(
            "neighbors_vs_cosine_similarity.png",
            bbox_extra_artists=(legend,),
            bbox_inches="tight",
        )
        plt.close()

        # save all of the files as csv with numpy
        import numpy as np

        np.savetxt(
            "minimal_cosine_similarity_in_circle.csv",
            np.array(minimal_cosine_similarity_in_circle),
            delimiter=",",
        )
        import pandas as pd

        # Create a combined DataFrame with all data
        combined_df = pd.DataFrame(
            {
                "index": range(len(number_of_neighbours)),
                "protein_id": [
                    list(data_working.values())[i]["protein_id"]
                    if i < len(data_working)
                    else None
                    for i in range(len(number_of_neighbours))
                ],
                "family_name": [
                    list(data_working.values())[i]["family_name"]
                    if i < len(data_working)
                    else None
                    for i in range(len(number_of_neighbours))
                ],
                "number_of_neighbours": [
                    ",".join(map(str, row))
                    if isinstance(row, (list, tuple))
                    else str(row)
                    for row in number_of_neighbours
                ],
                "cosine_similarity_scores": [
                    ",".join(map(str, row))
                    if isinstance(row, (list, tuple))
                    else str(row)
                    for row in cosine_similarity_scores_vs_neighbors_data_lines
                ],
                "ids_in_circle": ids_in_circle,
                "minimal_cosine_similarity": minimal_cosine_similarity_in_circle,
            }
        )

        # Save the combined DataFrame
        combined_df.to_csv("combined_protein_data.csv", index=False)

        np.savetxt(
            "already_processed_ids.csv",
            np.array(already_processed_ids),
            delimiter=",",
            fmt="%s",
        )

    # nohup python scr/code/013_TEM_Circular_Detection.py > 013_TEM_Circular_Detection.log 2>&1 &
