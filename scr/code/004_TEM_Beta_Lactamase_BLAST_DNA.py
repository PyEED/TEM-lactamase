# this is the big data pull for the beta lactamase data
# everything works now blast and backup so big stable data pull

# ------------------------------------- IMPORTING PACKAGES -------------------------------------
import logging
import os
from datetime import datetime

from dotenv import load_dotenv
from pyeed import Pyeed
from pyeed.tools.blast import Blast

# ------------------------------------- SETUP -------------------------------------

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
LOGGER = logging.getLogger(__name__)


uri = "bolt://129.69.129.130:7687"
user = "neo4j"
eedb = Pyeed(uri, user=user, password='12345678')
eedb.db.initialize_db_constraints(user, '12345678')

path_to_data = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull"
path_to_data_blast_dna = (
    "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/blast_data_dna"
)
path_to_data_backup = "/home/nab/Niklas/TEM-lactamase/data/003_data_pull/backup_data"
path_to_db_blast = "/databases/nt"

blast = Blast(
    mode="blastn",
    db_path=path_to_db_blast,
    db_name="nt",
    evalue=0.001,
    max_target_seqs=5000,
)

# ------------------------------------- FUNCTIONS -------------------------------------


def run_blast_and_fetch_data_dnas(id, path_to_data_blast_dna):
    """
    This function is supposed to identify all the dna sequences in the pull proteins in the database
    These DNA sequences are then blasted against the nt database. We want those blast results and save them to a csv file.
    """
    eedb.fetch_from_primary_db(id, db="ncbi_nucleotide")


    query_sequence = f"""
        MATCH (d:DNA) WHERE d.accession_id = "{id}"
        RETURN d.sequence
    """

    sequence = eedb.db.execute_read(query_sequence)[0]['d.sequence']

    df_blast = blast.search(sequence)

    print(df_blast.head())
    df_blast.to_csv(f"{path_to_data_blast_dna}/{id}.csv", index=False)


if __name__ == "__main__":

    # create folder of the timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    path_to_data_blast_dna = f"{path_to_data_blast_dna}/{timestamp}"
    path_to_data_backup = f"{path_to_data_backup}/{timestamp}"

    os.makedirs(path_to_data_blast_dna, exist_ok=True)
    os.makedirs(path_to_data_backup, exist_ok=True)

    query_ids_of_dna_connected_to_proteins = """
        MATCH (d:DNA)-[:ENCODES]->(p:Protein) RETURN d.accession_id
    """

    l = ['DAGOQK010000132.1', 'KT400838.1', 'KY727696.1', 'KT408379.1', 'KY719122.1', 'KY716156.1', 'KY718613.1', 'MG674491.1', 'KY713825.1', 'KY716608.1', 'KY717934.1', 'KY717430.1', 'KY727945.1', 'KY716813.1', 'KY733260.1', 'KY714344.1', 'KY733964.1', 'KY716804.1', 'KY716794.1', 'KY733612.1', 'KY714315.1', 'DAJIMF010000130.1', 'KY716589.1', 'KY716563.1', 'KY717862.1', 'KY713630.1', 'KY722742.1', 'ABFWUC010000143.1', 'KY716727.1', 'KY723904.1', 'KY719347.1', 'KY732751.1', 'DABEPA010000111.1', 'KY720031.1', 'KY717620.1', 'KY716595.1', 'KY717208.1', 'KY718331.1', 'KY718348.1', 'KT416718.1', 'KY718519.1', 'KY713818.1', 'KY713793.1', 'KY713894.1', 'KY724709.1', 'KY717113.1', 'KY716564.1', 'KY713625.1', 'KY716978.1', 'KY713638.1', 'KY713730.1', 'KY718112.1', 'KM017939.1', 'KY717286.1', 'DAFUDQ010000307.1', 'KY731731.1', 'KY716593.1', 'KT403343.1', 'KY718632.1', 'KY723927.1', 'KU664682.1', 'KY717099.1', 'KT407978.1', 'ABFXVI020000093.1', 'KY717922.1', 'KY716580.1', 'KY733303.1', 'KT400469.1', 'KT417086.1', 'KY714122.1', 'KT403629.1', 'ABLWLE010000045.1', 'KY717499.1', 'KY718113.1', 'KY716761.1', 'KY727935.1', 'KY713650.1', 'KY724696.1', 'KY732641.1', 'KY717285.1', 'KY732132.1', 'KY713861.1', 'KY728071.1', 'KT395322.1', 'KY717047.1', 'KY719110.1', 'KY715784.1', 'KY716929.1', 'KY727641.1', 'KY718292.1', 'KY717226.1', 'KY714225.1', 'KY739647.1', 'KY717543.1', 'DAGORP010000051.1', 'KY713847.1', 'KY717739.1', 'KY717390.1', 'KY724298.1', 'KY713759.1', 'KY727652.1', 'KY717719.1', 'KY716710.1', 'KY716971.1', 'KY716543.1', 'ABFBHS010000207.1', 'KY715060.1', 'KY718557.1', 'ABGJUB010000054.1', 'KY719130.1', 'KY713651.1', 'KT395241.1', 'KY732403.1']

    query_ids_of_dna_connected_to_proteins = eedb.db.execute_read(
        query_ids_of_dna_connected_to_proteins
    )

    print(
        f"Number of DNA sequences to blast: {len(query_ids_of_dna_connected_to_proteins)- len(l)}"
    )

    for id in query_ids_of_dna_connected_to_proteins:
        if id["d.accession_id"] in l:
            continue
        id_clean = id["d.accession_id"]
        print(f"Blasting {id_clean}")
        run_blast_and_fetch_data_dnas(id_clean, path_to_data_blast_dna)
        break


# nohup python scr/code/004_TEM_Beta_Lactamase_BLAST_DNA.py > 004_TEM_Beta_Lactamase_BLAST_DNA.log 2>&1 &
