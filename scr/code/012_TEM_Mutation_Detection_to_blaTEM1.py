import logging
import os

from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.embedding_analysis import EmbeddingTool
from pyeed.analysis.mutation_detection import MutationDetection
from pyeed.analysis.standard_numbering import StandardNumberingTool

path_to_data_blast = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna/2025-01-19_12-37-48"


load_dotenv()
password = os.getenv("NEO4J_NIKLAS_TEM_NEW_START")
if password is None:
    raise ValueError("KEY is not set in the .env file.")


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)


uri = "bolt://129.69.129.130:2127"
user = "neo4j"
eedb = Pyeed(uri, user=user, password=password)
eedb.db.initialize_db_constraints(user, password)


# we start by reading in all 258 TEM-lactamase proteins and check their identical ids
ids_tem = {}

base_url_tem_family_card = "http://purl.obolibrary.org/obo/ARO_3000014"

# get all the children of the TEM-lactamase family
query = f"""
MATCH (o:OntologyObject {{name: '{base_url_tem_family_card}'}})-[*1..1]-(n) RETURN n
"""

result = eedb.db.execute_read(query)

for single_tem in result:
    if single_tem["n"]["name"] == "http://purl.obolibrary.org/obo/ARO_3000078":
        continue

    tem_name = single_tem["n"]["label"]
    tem_url = single_tem["n"]["name"]
    ids_tem[tem_name] = {"tem_name": tem_name, "tem_url": tem_url}

    # now we check for the URL and get the matching protein and read out the number of IdenticalIds
    query_tem_url = f"""
    MATCH (o:OntologyObject {{name: '{tem_url}'}})-[*1..1]-(n:Protein) RETURN n
    """

    result_tem_url = eedb.db.execute_read(query_tem_url)
    if len(result_tem_url) == 0:
        continue
    result_tem_url = result_tem_url[0]

    if "IdenticalIds" in result_tem_url["n"]:
        ids_tem[tem_name]["identical_ids"] = result_tem_url["n"]["IdenticalIds"]
    else:
        ids_tem[tem_name]["identical_ids"] = []

    ids_tem[tem_name]["accession_id"] = result_tem_url["n"]["accession_id"]

print(f"The number of TEM-lactamase proteins is: {len(ids_tem)}")
print(ids_tem)

# max number of neighbours
n_neighbours = 500
name_of_standard_numbering_tool = (
    "standard_numbering_pairwise_blaTEM1a_mutations_to_blaTEM1"
)

et = EmbeddingTool()
sn = StandardNumberingTool(name=name_of_standard_numbering_tool)
md = MutationDetection()

blaTEM1a_id = "CAD09800.1"
blaTEM1a_database_id = None


# find the coresponding database id
for tem_name, tem_data in ids_tem.items():
    if "accession_id" in tem_data:
        if tem_data["accession_id"] == blaTEM1a_id:
            blaTEM1a_database_id = tem_data["accession_id"]
            break
        else:
            if blaTEM1a_id in tem_data["identical_ids"]:
                blaTEM1a_database_id = tem_data["accession_id"]
                break

print(blaTEM1a_database_id)


if __name__ == "__main__":
    already_processed_pairs = []

    for tem_name, tem_data in ids_tem.items():
        if "accession_id" in tem_data:
            # get the closest neighbours
            results = et.find_nearest_neighbors_based_on_vector_index(
                index_name="vector_index_Protein_embedding",
                query_protein_id=tem_data["accession_id"],
                number_of_neighbors=n_neighbours,
                db=eedb.db,
            )

            ids = [neighbour[0] for neighbour in results] + [tem_data["accession_id"]]

            sn.apply_standard_numbering_pairwise(
                base_sequence_id=blaTEM1a_database_id, db=eedb.db, list_of_seq_ids=ids
            )

            # we need to create all of the permutations of the neighbours with the base sequence
            # please that the reverse direction should not be included
            # this means that the base sequence is always the first element in the tuple and the second element is the neighbour
            permutations = [(blaTEM1a_database_id, neighbour) for neighbour in ids]
            # print(f"The permutations of the neighbours including the base sequence are: {len(permutations)}")

            # we now want to exclude the pairs that we already processed keeping in mind that we always add in the list both directions
            permuations_to_process = [
                pair for pair in permutations if pair not in already_processed_pairs
            ]
            LOGGER.info(
                f"The number of permutations to process is: {len(permuations_to_process)}"
            )

            # we now update the already_processed_pairs list with the new pairs
            # we need to add the reverse of the pair as well
            already_processed_pairs.extend(
                [(pair[1], pair[0]) for pair in permuations_to_process]
            )
            already_processed_pairs.extend(permuations_to_process)

            for pair in permuations_to_process:
                if pair[0] == pair[1]:
                    continue

                LOGGER.info(f"Processing pair {pair[0]} and {pair[1]}")

                mutations = md.get_mutations_between_sequences(
                    pair[0],
                    pair[1],
                    eedb.db,
                    name_of_standard_numbering_tool,
                    save_to_db=True,
                )

                LOGGER.info(
                    f"The mutations are: {mutations}, there are {len(mutations['from_positions'])} mutations"
                )

    # nohup python scr/code/012_TEM_Mutation_Detection_to_blaTEM1.py > 012_TEM_Mutation_Detection_to_blaTEM1.log 2>&1 &
