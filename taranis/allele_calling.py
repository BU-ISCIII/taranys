import io
import logging
import os
import rich.console


import taranis.utils
import taranis.blast
import pdb

# import numpy
from collections import OrderedDict
from pathlib import Path
from Bio import SeqIO


log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)


class AlleleCalling:
    def __init__(
        self, sample_file: str, schema: str, reference_alleles: list, out_folder: str
    ):
        # self.prediction = prediction
        self.sample_file = sample_file
        self.schema = schema
        self.ref_alleles = reference_alleles
        self.out_folder = out_folder
        self.s_name = Path(sample_file).stem
        self.blast_dir = os.path.join(out_folder, "blastdb")
        # create blast for sample file
        self.blast_obj = taranis.blast.Blast("nucl")
        _ = self.blast_obj.create_blastdb(sample_file, self.blast_dir)

    def assign_allele_type(
        self, blast_result: list, allele_file: str, allele_name: str
    ) -> list:
        def get_blast_details(blast_result: list, allele_name: str) -> list:
            # get blast details
            blast_details = [
                blast_result[0].split("_")[0],  # Core gene
                self.s_name,
                "gene annotation",
                "product annotation",
                allele_name,
                "allele quality",
                blast_result[1],  # contig
                blast_result[3],  # query length
                blast_result[9],  # contig start
                blast_result[10],  # contig end
                blast_result[13],  # contig sequence
            ]
            return blast_details

        if len(blast_result) > 1:
            # allele is named as NIPHEM
            column_blast_res = blast_result[0].split("\t")
            column_blast_res[13] = column_blast_res[13].replace("-", "")
            allele_details = get_blast_details(column_blast_res, allele_name)
            return ["NIPHEM", allele_name, allele_details]


        elif len(blast_result) == 1:
            column_blast_res = blast_result[0].split("\t")
            column_blast_res[13] = column_blast_res[13].replace("-", "")
            allele_details = get_blast_details(column_blast_res, allele_name)

            grep_result = taranis.utils.grep_execution(
                allele_file, column_blast_res[13], "-b1"
            )
            # check if sequence match alleles in schema
            if len(grep_result) > 0:
                allele_name = grep_result[0].split(">")[1]

                # allele is labled as EXACT
                return ["EXC", allele_name, allele_details]
            # check if contig is shorter than allele
            if int(column_blast_res[3]) > int(column_blast_res[4]):
                # check if sequence is shorter because it starts or ends at the contig
                if (
                    column_blast_res[9] == 1  # check  at contig start
                    or column_blast_res[14]
                    == column_blast_res[10]  # check at contig end
                    or column_blast_res[10] == 1  # check reverse at contig end
                    or column_blast_res[9]
                    == column_blast_res[14]  # check reverse at contig start
                ):
                    # allele is labled as PLOT
                    return ["PLOT", allele_name, allele_details]
                # allele is labled as ASM
                return ["ASM", allele_name, allele_details]
            # check if contig is longer than allele
            if int(column_blast_res[3]) < int(column_blast_res[4]):
                # allele is labled as ALM
                return ["ALM", allele_name, allele_details]
            if int(column_blast_res[3]) == int(column_blast_res[4]):
                # allele is labled as INF
                return ["INF", allele_name, allele_details]
        else:
            print("ERROR: No blast result found")
            return ["LNF", "allele_name", "LNF"]

    def search_match_allele(self):
        # Create  blast db with sample file

        result = {"allele_type": {}, "allele_match": {}, "allele_details": {}}
        count = 0
        for ref_allele in self.ref_alleles:
            count += 1
            print(" Processing allele ", ref_allele, " ", count, " of ", len(self.ref_alleles))
            # schema_alleles = os.path.join(self.schema, ref_allele)
            # parallel in all CPUs in cluster node
            alleles = OrderedDict()
            match_found = False
            with open(ref_allele, "r") as fh:
                for record in SeqIO.parse(fh, "fasta"):
                    alleles[record.id] = str(record.seq)
            count_2 = 0
            for r_id, r_seq in alleles.items():
                count_2 += 1

                print("Running blast for ", count_2, " of ", len(alleles))
                # create file in memory to increase speed
                query_file = io.StringIO()
                query_file.write(">" + r_id + "\n" + r_seq)
                query_file.seek(0)
                blast_result = self.blast_obj.run_blast(
                    query_file.read(), perc_identity=90, num_threads=4, query_type="stdin"
                )
                if len(blast_result) > 0:
                    match_found = True
                    break
            # Close object and discard memory buffer
            query_file.close()
            if match_found:
                allele_file = os.path.join(self.schema, os.path.basename(ref_allele))
                # blast_result = self.blast_obj.run_blast(q_file,perc_identity=100)
                try:
                    allele_name = Path(allele_file).stem
                    (
                        result["allele_type"][allele_name],
                        result["allele_match"][allele_name],
                        result["allele_details"][allele_name],
                    ) = self.assign_allele_type(blast_result, allele_file, allele_name)
                except Exception as e:
                    stderr.print(f"Error: {e}")

            else:
                # Sample does not have a reference allele to be matched
                # Keep LNF info
                # ver el codigo de espe
                # lnf_tpr_tag()
                pass

        return result


def parallel_execution(
    sample_file: str, schema: str, reference_alleles: list, out_folder: str
):
    allele_obj = AlleleCalling(sample_file, schema, reference_alleles, out_folder)
    return allele_obj.search_match_allele()
