# from Bio.Seq import Seq

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess

# import taranys.utils

import random

"""
    Para hacer las pruebas con alfaclust activo el entorno de conda alfatclust_env
    despues me voy a la carpeta donde me he descargado, de git, alfatclust y
    ejecuto :
     ./alfatclust.py -i /media/lchapado/Reference_data/proyectos_isciii/taranys/taranys_testing_data/listeria_testing_schema/lmo0003.fasta  -o /media/lchapado/Reference_data/proyectos_isciii/taranys/test/alfacluster_test/resultado_alfaclust_lmo003  -l 0.9
    despues ejecuto este programa de prueba cambiando los ficheros de resultados

"""

# read result of alfatclust

alfa_clust_file = "/media/lchapado/Reference_data/proyectos_isciii/taranys/test/resultado_alfatclust-090"
with open(alfa_clust_file, "r") as fh:
    lines = fh.readlines()
alleles_found = False
locus_list = []
for line in lines:
    line = line.strip()
    if line == "#Cluster 5":
        if alleles_found is False:
            alleles_found = True
            continue
    if alleles_found:
        if "#Cluster" in line:
            break
        locus_list.append(line)

rand_locus = random.choice(locus_list)
schema_file = "/media/lchapado/Reference_data/proyectos_isciii/taranys/taranys_testing_data/listeria_testing_schema/lmo0002.fasta"
new_schema_file = (
    "/media/lchapado/Reference_data/proyectos_isciii/taranys/test/cluster_lmo0002.fasta"
)
q_file = "/media/lchapado/Reference_data/proyectos_isciii/taranys/test/q_file.fasta"
with open(schema_file) as fh:
    with open(new_schema_file, "w") as fo:
        for record in SeqIO.parse(schema_file, "fasta"):
            if record.id in locus_list:
                SeqIO.write(record, fo, "fasta")

# choose a random locus for testing
with open(new_schema_file) as fh:
    with open(q_file, "w") as fo:
        for record in SeqIO.parse(new_schema_file, "fasta"):
            if record.id == rand_locus:
                SeqIO.write(record, fo, "fasta")
                break
print("Selected locus: ", rand_locus)
db_name = "/media/lchapado/Reference_data/proyectos_isciii/taranys/test/testing_clster/lmo0002"
blast_command = [
    "makeblastdb",
    "-in",
    new_schema_file,
    "-parse_seqids",
    "-dbtype",
    "nucl",
    "-out",
    db_name,
]
blast_result = subprocess.run(
    blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
)

blast_parameters = '"6 , qseqid , sseqid , pident ,  qlen , length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend , sseq , qseq"'
# db=self.blast_dir, evalue=evalue, perc_identity=perc_identity_ref, reward=reward, penalty=penalty, gapopen=gapopen, gapextend=gapextend, outfmt=blast_parameters, max_target_seqs=max_target_seqs, max_hsps=max_hsps, num_threads=num_threads, query=core_reference_allele_path)
cline = NcbiblastnCommandline(
    db=db_name,
    evalue=0.001,
    perc_identity=90,
    reward=1,
    penalty=-2,
    gapopen=1,
    gapextend=1,
    outfmt=blast_parameters,
    max_target_seqs=1100,
    max_hsps=1000,
    num_threads=4,
    query=q_file,
)

try:
    out, _ = cline()
except Exception as e:
    print(e)
b_lines = out.splitlines()
print("longitud del cluster = ", len(locus_list))
print("numero de matches = ", len(b_lines))
