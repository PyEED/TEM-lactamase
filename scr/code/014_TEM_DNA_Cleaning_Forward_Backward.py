import logging
import os

from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.analysis.embedding_analysis import EmbeddingTool
from pyeed.analysis.mutation_detection import MutationDetection
from pyeed.analysis.sequence_alignment import PairwiseAligner

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


et = EmbeddingTool()
md = MutationDetection()
pa = PairwiseAligner(node_type="DNA")


blaTEM1a_id = "AAB59737.1"  # Lohey
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


def complement_dna_sequence(dna_sequence):
    # complement the dna sequence
    dna_sequence = list(dna_sequence)
    dna_sequence.reverse()
    dna_sequence = "".join(dna_sequence)
    return (
        dna_sequence.replace("A", "t")
        .replace("T", "a")
        .replace("C", "g")
        .replace("G", "c")
        .upper()
    )


def count_number_of_completed_regions(accession_id_blaTEM1a):
    query_count_number_of_completed_regions = """
    MATCH (r1:Region)-[r1d1:HAS_REGION]-(d1:DNA)-[:ENCODES]-(p1:Protein)-[:MUTATION]-(p2:Protein)-[:ENCODES]-(d2:DNA)-[r2d2:HAS_REGION]-(r2:Region)
    WHERE p1.accession_id = $accession_id_blaTEM1a
      AND r1.sequence_id = p1.accession_id
      AND d2.complement = 1
    RETURN count(DISTINCT d2)
    """

    result_count_number_of_completed_regions = eedb.db.execute_read(
        query_count_number_of_completed_regions,
        parameters={"accession_id_blaTEM1a": accession_id_blaTEM1a},
    )

    return result_count_number_of_completed_regions[0]["count(DISTINCT d2)"]


def make_reset_of_completed_regions(accession_id_blaTEM1a):
    query_make_reset_of_completed_regions = """
    MATCH (r1:Region)-[r1d1:HAS_REGION]-(d1:DNA)-[:ENCODES]-(p1:Protein)-[:MUTATION]-(p2:Protein)-[:ENCODES]-(d2:DNA)-[r2d2:HAS_REGION]-(r2:Region)
    WHERE p1.accession_id = $accession_id_blaTEM1a
      AND r1.sequence_id = p1.accession_id
      AND d2.completed = 1
    WITH DISTINCT d2, r2d2, r2d2.old_start as original_start, r2d2.old_end as original_end
    SET d2.sequence = d2.old_sequence,
        r2d2.start = original_start,
        r2d2.end = original_end
    REMOVE d2.completed, d2.old_sequence, d2.complement,
           r2d2.old_start, r2d2.old_end
    """


if __name__ == "__main__":
    key_base_sequence_id = "AL513383.1"
    key_base_sequence_Region_id = 10228113

    print(
        count_number_of_completed_regions(
            accession_id_blaTEM1a=blaTEM1a_database_id,
        )
    )

    # find all Protein and their DNA which are not yet connected through a region region mutation to the base sequence DNA
    # the idea is that the protein are already connected and now I want to find the dna to the proteins whoch are already connected and get their region ids

    # The previous query might return duplicate entries because:
    # 1. Multiple paths between proteins
    # 2. Multiple regions for the same DNA
    # Using DISTINCT to ensure unique region and DNA combinations
    query_of_standard_numbering_tool_dna_regions_ids = """
    MATCH (r1:Region)-[r1d1:HAS_REGION]-(d1:DNA)-[:ENCODES]-(p1:Protein)-[:MUTATION]-(p2:Protein)-[:ENCODES]-(d2:DNA)-[r2d2:HAS_REGION]-(r2:Region)
    WHERE p1.accession_id = $accession_id_blaTEM1a 
        AND r2.sequence_id = p2.accession_id 
        AND r1.annotation = $region_annotation
        AND r1.sequence_id = $accession_id_blaTEM1a
        AND r2.annotation = $region_annotation
        AND d1.accession_id = $key_base_sequence_id
        AND NOT (r1)-[:MUTATION]-(r2)
    RETURN id(r2), d2.accession_id, r1d1.start, r1d1.end, r2d2.start, r2d2.end, d1.sequence, d2.sequence, d1.accession_id
    """

    result_not_yet_mutated_regions = eedb.db.execute_read(
        query_of_standard_numbering_tool_dna_regions_ids,
        parameters={
            "accession_id_blaTEM1a": blaTEM1a_database_id,
            "region_annotation": "coding sequence",
            "key_base_sequence_id": key_base_sequence_id,
        },
    )

    print(len(result_not_yet_mutated_regions))

    region_ids = [
        result_not_yet_mutated_regions[i]["id(r2)"]
        for i in range(len(result_not_yet_mutated_regions))
    ]
    dna_accession_ids = [
        result_not_yet_mutated_regions[i]["d2.accession_id"]
        for i in range(len(result_not_yet_mutated_regions))
    ]

    region_ids_1_start_end = [
        (
            region_ids[i],
            result_not_yet_mutated_regions[i]["r1d1.start"],
            result_not_yet_mutated_regions[i]["r1d1.end"],
        )
        for i in range(len(region_ids))
    ]

    region_ids_2_start_end = [
        (
            region_ids[i],
            result_not_yet_mutated_regions[i]["r2d2.start"],
            result_not_yet_mutated_regions[i]["r2d2.end"],
        )
        for i in range(len(region_ids))
    ]

    dna_sequences = [
        (
            result_not_yet_mutated_regions[i]["d1.sequence"],
            result_not_yet_mutated_regions[i]["d2.sequence"],
        )
        for i in range(len(region_ids))
    ]

    # between those DNAs and the spefic regiosn referenced by the region i want a pairwise alignement done
    # if the pairwise alignemnet has a identity of lower than 95% i want to do it also backwards the dna complemented
    # then resulting from the two alignement the one with the higher identity should be written in the DNA at the sequence attribute the old sequence will move to the DNA attribute old_Sequence and a flag of completed should be set to 1
    # the pairwise aligner should be used from the pyeed package

    for i in range(len(region_ids)):
        alignment = pa.align_pairwise(
            seq1={
                key_base_sequence_id: dna_sequences[i][0][
                    region_ids_1_start_end[i][1] : region_ids_1_start_end[i][2]
                ],
            },
            seq2={
                dna_accession_ids[i]: dna_sequences[i][1][
                    region_ids_2_start_end[i][1] : region_ids_2_start_end[i][2]
                ],
            },
        )

        if alignment["identity"] < 0.95:
            alignment_complement = pa.align_pairwise(
                seq1={
                    key_base_sequence_id: dna_sequences[i][0][
                        region_ids_1_start_end[i][1] : region_ids_1_start_end[i][2]
                    ],
                },
                seq2={
                    dna_accession_ids[i]: complement_dna_sequence(
                        dna_sequences[i][1][
                            region_ids_2_start_end[i][1] : region_ids_2_start_end[i][2]
                        ]
                    ),
                },
            )

            if alignment_complement["identity"] > alignment["identity"]:
                # update the dna sequence
                print(len(dna_sequences[i][1]))
                dna_sequence_new = complement_dna_sequence(dna_sequences[i][1])
                dna_sequence_old = dna_sequences[i][1]

                old_start = region_ids_2_start_end[i][1]
                old_end = region_ids_2_start_end[i][2]

                sequence_length = len(dna_sequences[i][1])
                new_start = sequence_length - old_end
                new_end = sequence_length - old_start

                # update the dna sequence in the database
                query_update_dna_sequence = """
                MATCH (r1:Region {annotation: $region_annotation})-[rel_reg:HAS_REGION]-(d:DNA)-[:ENCODES]-(p:Protein)
                WHERE d.accession_id = $accession_id_dna AND id(r1) = $region_id
                SET d.sequence = $dna_sequence_new, d.old_sequence = $dna_sequence_old, d.complement = 1, rel_reg.start = $new_start, rel_reg.end = $new_end, rel_reg.old_start = $old_start, rel_reg.old_end = $old_end
                """

                eedb.db.execute_write(
                    query_update_dna_sequence,
                    parameters={
                        "accession_id_dna": dna_accession_ids[i],
                        "dna_sequence_new": dna_sequence_new,
                        "dna_sequence_old": dna_sequence_old,
                        "region_id": region_ids[i],
                        "region_annotation": "coding sequence",
                        "new_start": new_start,
                        "new_end": new_end,
                        "old_start": old_start,
                        "old_end": old_end,
                    },
                )

                print(
                    f"\n \n \n \n \n \n Updated DNA sequence for {dna_accession_ids[i]} with identity complement {alignment_complement['identity']} and original {alignment['identity']} and \n \n sequence new: {dna_sequence_new[new_start:new_end]} \n \n  sequence old: {dna_sequence_old[old_start:old_end]} \n \n sequence target: {dna_sequences[i][0][region_ids_1_start_end[i][1] : region_ids_1_start_end[i][2]]}"
                )


# nohup python scr/code/014_TEM_DNA_Cleaning_Forward_Backward.py > 014_TEM_DNA_Cleaning_Forward_Backward.log 2>&1 &