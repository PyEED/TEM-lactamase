import logging
import os

from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.embedding_analysis import EmbeddingTool
from pyeed.analysis.mutation_detection import MutationDetection
from pyeed.analysis.standard_numbering import StandardNumberingTool

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


# max number of neighbours
n_neighbours = 100

et = EmbeddingTool()
md = MutationDetection()

blaTEM1a_id = "WP_000027057.1"  # Lahey
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
    key_base_sequence_id = "AL513383.1"
    key_base_sequence_Region_id = 10228113

    name_of_standard_numbering_tool = (
        "standard_numbering_pairwise_blaTEM1a_DNA_region_10228113"
    )

    sn_dna = StandardNumberingTool(name=name_of_standard_numbering_tool)

    query_of_standard_numbering_tool_dna_regions_ids = """
    MATCH (p:Protein {accession_id: $accession_id_blaTEM1a})-[*1..1]-(d:DNA)-[rel:HAS_REGION]->(r:Region)
    WHERE r.annotation = $region_annotation AND r.sequence_id = $accession_id_blaTEM1a
    RETURN id(r), d.accession_id
    """

    result_dna_accession_ids = eedb.db.execute_read(
        query_of_standard_numbering_tool_dna_regions_ids,
        parameters={
            "accession_id_blaTEM1a": blaTEM1a_database_id,
            "region_annotation": "coding sequence",
        },
    )

    result_dna_accession_ids_accession_ids = [
        result_dna_accession_ids[i]["d.accession_id"]
        for i in range(len(result_dna_accession_ids))
    ]

    result_region_ids = [
        result_dna_accession_ids[i]["id(r)"]
        for i in range(len(result_dna_accession_ids))
    ]

    print(result_dna_accession_ids_accession_ids)
    print(result_region_ids)

    # apply the standard numbering tool to the DNAs
    sn_dna.apply_standard_numbering_pairwise(
        base_sequence_id=key_base_sequence_id,
        db=eedb.db,
        list_of_seq_ids=result_dna_accession_ids_accession_ids,
        node_type="DNA",
        region_ids_neo4j=result_region_ids,
    )

    # apply the mutation detection tool to the DNAs
    pairs = [(key_base_sequence_id, i) for i in result_dna_accession_ids_accession_ids]
    for pair in pairs:
        md.get_mutations_between_sequences(
            sequence_id1=pair[0],
            sequence_id2=pair[1],
            db=eedb.db,
            node_type="DNA",
            standard_numbering_tool_name=name_of_standard_numbering_tool,
            region_ids_neo4j=result_region_ids,
        )

    # find all Protein and their DNA which are not yet connected through a region region mutation to the base sequence DNA
    # the idea is that the protein are already connected and now I want to find the dna to the proteins whoch are already connected and get their region ids

    # The previous query might return duplicate entries because:
    # 1. Multiple paths between proteins
    # 2. Multiple regions for the same DNA
    # Using DISTINCT to ensure unique region and DNA combinations
    query_of_standard_numbering_tool_dna_regions_ids = """
    MATCH (r1:Region)-[:HAS_REGION]-(d1:DNA)-[:ENCODES]-(p1:Protein)-[:MUTATION]-(p2:Protein)-[:ENCODES]-(d2:DNA)-[:HAS_REGION]-(r2:Region)
    WHERE p1.accession_id = $accession_id_blaTEM1a 
      AND r2.sequence_id = p2.accession_id 
      AND r1.annotation = $region_annotation 
      AND r2.annotation = $region_annotation
      AND NOT (r1)-[:MUTATION]-(r2)
    RETURN DISTINCT id(r2), d2.accession_id
    """

    result_not_yet_mutated_regions = eedb.db.execute_read(
        query_of_standard_numbering_tool_dna_regions_ids,
        parameters={
            "accession_id_blaTEM1a": blaTEM1a_database_id,
            "region_annotation": "coding sequence",
        },
    )

    region_ids = [
        result_not_yet_mutated_regions[i]["id(r2)"]
        for i in range(len(result_not_yet_mutated_regions))
    ]
    dna_accession_ids = [
        result_not_yet_mutated_regions[i]["d2.accession_id"]
        for i in range(len(result_not_yet_mutated_regions))
    ]

    print(len(region_ids))  # 4801

    batch_size = 10
    # in batches of 1000 apply these region ids to the standard numbering tool and detect the mutations
    for i in range(0, len(region_ids), batch_size):
        sn_dna.apply_standard_numbering_pairwise(
            base_sequence_id=key_base_sequence_id,
            db=eedb.db,
            list_of_seq_ids=dna_accession_ids[i : i + batch_size],
            node_type="DNA",
            region_ids_neo4j=region_ids[i : i + batch_size]
            + [key_base_sequence_Region_id],
        )

        # apply the mutation detection tool to the DNAs
        for j, accession_id in enumerate(dna_accession_ids[i : i + batch_size]):
            md.get_mutations_between_sequences(
                sequence_id1=key_base_sequence_id,
                sequence_id2=accession_id,
                db=eedb.db,
                node_type="DNA",
                standard_numbering_tool_name=name_of_standard_numbering_tool,
                region_ids_neo4j=region_ids[i + j : i + j + 1]
                + [key_base_sequence_Region_id],
            )
