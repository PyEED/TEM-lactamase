import logging
import os

from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.embedding_analysis import EmbeddingTool
from pyeed.analysis.mutation_detection import MutationDetection
from pyeed.analysis.standard_numbering import StandardNumberingTool

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


# max number of neighbours
n_neighbours = 100
name_of_standard_numbering_tool = "standard_numbering_pairwise_blaTEM1a_DNA_KY739721.1"

et = EmbeddingTool()
sn_dna = StandardNumberingTool(name=name_of_standard_numbering_tool)
md = MutationDetection()

blaTEM1a_id = "AAB59737.1"
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
    # apply the standard numbering tool to the DNAs that are connected to the blaTEM1a protein
    query = f"""
    MATCH (p:Protein {{accession_id: '{blaTEM1a_database_id}'}})-[*1..1]-(d:DNA) RETURN d
    """

    result = eedb.db.execute_read(query)

    result_dna_accession_ids = []
    for single_dna in result:
        result_dna_accession_ids.append(single_dna["d"]["accession_id"])

    print(result_dna_accession_ids)

    # apply the standard numbering tool to the DNAs
    sn_dna.apply_standard_numbering_pairwise(
        base_sequence_id="KY739721.1",
        db=eedb.db,
        list_of_seq_ids=result_dna_accession_ids,
        node_type="DNA",
    )

    # apply the mutation detection tool to the DNAs
    pairs = [("KY739721.1", i) for i in result_dna_accession_ids]
    for pair in pairs:
        md.get_mutations_between_sequences(
            sequence_id1=pair[0],
            sequence_id2=pair[1],
            db=eedb.db,
            node_type="DNA",
            standard_numbering_tool_name=name_of_standard_numbering_tool,
        )
