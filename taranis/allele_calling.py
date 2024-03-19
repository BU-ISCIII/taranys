import io
import logging
import os
import rich.console

import taranis.utils
import taranis.blast
import pdb
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
        self,
        sample_file: str,
        schema: str,
        annotation: dict,
        reference_alleles: list,
        out_folder: str,
        inf_alle_obj: object,
        snp_request: bool = False,
        aligment_request: bool = False,
    ):
        self.prediction_data = annotation  # store prediction annotation
        self.sample_file = sample_file
        self.schema = schema
        self.ref_alleles = reference_alleles
        self.out_folder = out_folder
        self.s_name = Path(sample_file).stem
        self.blast_dir = os.path.join(out_folder, "blastdb")
        # create blast for sample file
        self.blast_obj = taranis.blast.Blast("nucl")
        _ = self.blast_obj.create_blastdb(sample_file, self.blast_dir)
        # store inferred allele object
        self.inf_alle_obj = inf_alle_obj
        self.snp_request = snp_request
        self.aligment_request = aligment_request

    def assign_allele_type(
        self, blast_result: list, allele_file: str, allele_name: str
    ) -> list:
        def get_blast_details(blast_result: list, allele_name: str) -> list:
            match_allele_name = blast_result[0]
            try:
                gene_annotation = self.prediction_data[match_allele_name]["gene"]
                product_annotation = self.prediction_data[match_allele_name]["product"]
                allele_quality = self.prediction_data[match_allele_name][
                    "allele_quality"
                ]
            except KeyError:
                gene_annotation = "Not found"
                product_annotation = "Not found"
                allele_quality = "Not found"
            if int(blast_result[10]) > int(blast_result[9]):
                direction = "+"
            else:
                direction = "-"
            # get blast details
            blast_details = [
                self.s_name,  # sample name
                blast_result[1],  # contig
                allele_name,  # core gene name
                blast_result[0],  # allele gene
                "coding",  # coding allele type. To be filled later idx = 4
                blast_result[3],  # query length
                blast_result[4],  # match length
                blast_result[14],  # contig length
                blast_result[9],  # contig start
                blast_result[10],  # contig end
                direction,
                gene_annotation,
                product_annotation,
                allele_quality,
                blast_result[13],  # contig sequence
            ]

            return blast_details

        valid_blast_result = []
        match_full_length = 0
        match_partial_length = 0
        for b_result in blast_result:
            column_blast_res = b_result.split("\t")
            query_length = int(column_blast_res[4]) / int(column_blast_res[3])
            if query_length >= 0.8:
                valid_blast_result.append(b_result)
                if query_length == 1:
                    match_full_length += 1
                else:
                    match_partial_length += 1

        if len(valid_blast_result) > 1:
            # allele could be named as NIPHEM or NIPH
            multi_allele = []
            for valid_result in valid_blast_result:
                column_blast_res = valid_result.split("\t")
                column_blast_res[13] = column_blast_res[13].replace("-", "")
                allele_details = get_blast_details(column_blast_res, allele_name)
                if match_full_length >= 2:
                    # labled as NIPHEM if all alleles are in the same contig
                    allele_details[4] = "NIPHEM_" + allele_details[3]
                    clasification = "NIPHEM"
                else:
                    # labled as NIPH if all alleles are in different contigs
                    allele_details[4] = "NIPH_" + allele_details[3]
                    clasification = "NIPH"
                multi_allele.append(allele_details)
            return [clasification, allele_name, multi_allele]

        elif len(valid_blast_result) == 1:
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
                allele_details[4] = "EXC_" + allele_details[3]
                return ["EXC", allele_name, allele_details]
            # check if contig is shorter than allele
            if int(column_blast_res[3]) > int(column_blast_res[4]):
                # check if sequence is shorter because it starts or ends at the contig
                if (
                    column_blast_res[9] == "1"  # check  at contig start
                    or column_blast_res[14]
                    == column_blast_res[10]  # check at contig end
                    or column_blast_res[10] == "1"  # check reverse at contig end
                    or column_blast_res[9]
                    == column_blast_res[14]  # check reverse at contig start
                ):
                    # allele is labled as PLOT
                    allele_details[4] = "PLOT_" + allele_details[3]
                    return ["PLOT", allele_name, allele_details]
                # allele is labled as ASM
                allele_details[4] = "ASM_" + allele_details[3]
                return ["ASM", allele_name, allele_details]
            # check if contig is longer than allele
            if int(column_blast_res[3]) < int(column_blast_res[4]):
                # allele is labled as ALM
                allele_details[4] = "ALM_" + allele_details[3]
                return ["ALM", allele_name, allele_details]
            if int(column_blast_res[3]) == int(column_blast_res[4]):
                # allele is labled as INF
                allele_details[4] = (
                    "INF_"
                    + allele_name
                    + "_"
                    + str(
                        self.inf_alle_obj.get_inferred_allele(
                            column_blast_res[14], allele_name
                        )
                    )
                )
                return ["INF", allele_name, allele_details]
        else:
            # analyze again the blast result to check with lower query size, 0.75
            # it starts/ends at the contig. Then it is labled as PLOT

            multi_allele = []
            clasification = ""
            for b_result in blast_result:
                column_blast_res = b_result.split("\t")
                query_length = int(column_blast_res[4]) / int(column_blast_res[3])
                if query_length >= 0.75:
                    if (
                        column_blast_res[9] == "1"  # check  at contig start
                        or column_blast_res[14]
                        == column_blast_res[10]  # check at contig end
                        or column_blast_res[10] == "1"  # check reverse at contig end
                        or column_blast_res[9]
                        == column_blast_res[14]  # check reverse at contig start
                    ):
                        allele_details = get_blast_details(
                            column_blast_res, allele_name
                        )
                        # allele is labled as PLOT
                        allele_details[4] = "PLOT_" + allele_details[3]
                        multi_allele.append(allele_details)
                        clasification = "PLOT"
            if clasification == "PLOT":
                return [clasification, allele_name, multi_allele]
            else:
                return ["LNF", "allele_name", "LNF"]

    def search_match_allele(self):
        # Create  blast db with sample file

        result = {
            "allele_type": {},
            "allele_match": {},
            "allele_details": {},
            "snp_data": {},
        }
        count = 0
        for ref_allele in self.ref_alleles:
            count += 1
            print(
                " Processing allele ",
                ref_allele,
                " ",
                count,
                " of ",
                len(self.ref_alleles),
            )
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
                    query_file.read(),
                    perc_identity=90,
                    num_threads=4,
                    query_type="stdin",
                )
                if len(blast_result) > 0:
                    match_found = True
                    break
            # Close object and discard memory buffer
            query_file.close()
            if match_found:
                allele_file = os.path.join(self.schema, os.path.basename(ref_allele))
                # blast_result = self.blast_obj.run_blast(q_file,perc_identity=100)
                allele_name = Path(allele_file).stem
                (
                    result["allele_type"][allele_name],
                    result["allele_match"][allele_name],
                    result["allele_details"][allele_name],
                ) = self.assign_allele_type(blast_result, allele_file, allele_name)
            else:
                # Sample does not have a reference allele to be matched
                # Keep LNF info
                result["allele_type"][allele_name] = "LNF"
                result["allele_match"][allele_name] = allele_name
                result["allele_details"][allele_name] = "LNF"
            if self.snp_request and result["allele_type"][allele_name] == "INF":
                # run snp analysis
                allele_seq = result["allele_details"][allele_name][14]
                result["snp_data"][allele_name] = taranis.utils.get_snp_position(
                    allele_seq, alleles
                )
        return result


def parallel_execution(
    sample_file: str,
    schema: str,
    prediction_data: dict,
    reference_alleles: list,
    out_folder: str,
    inf_alle_obj: object,
    snp_request: bool = False,
    aligment_request: bool = False,
):
    allele_obj = AlleleCalling(
        sample_file,
        schema,
        prediction_data,
        reference_alleles,
        out_folder,
        inf_alle_obj,
        snp_request,
        aligment_request,
    )
    return allele_obj.search_match_allele()


def collect_data(
    results: list, output: str, snp_request: bool, aligment_request: bool
) -> None:
    summary_result_file = os.path.join(output, "allele_calling_summary.csv")
    sample_allele_match_file = os.path.join(output, "allele_calling_match.csv")
    sample_allele_detail_file = os.path.join(output, "matching_contig.csv")
    allele_types = ["NIPHEM", "NIPH", "EXC", "PLOT", "ASM", "ALM", "INF", "LNF"]
    detail_heading = [
        "sample",
        "contig",
        "core gene",
        "allele name",
        "codification",
        "query lenght",
        "match lengt",
        "contig length",
        "contig start",
        "contig stop",
        "direction",
        "gene notation",
        "product notation",
        "allele quality",
        "sequence",
    ]
    summary_result = {}
    sample_allele_match = {}

    # get allele list
    first_sample = list(results[0].keys())[0]
    allele_list = sorted(results[0][first_sample]["allele_type"].keys())
    for result in results:
        for sample, values in result.items():
            sum_allele_type = OrderedDict()  # used for summary file
            allele_match = {}
            for allele_type in allele_types:
                sum_allele_type[allele_type] = 0
            for allele, type_of_allele in values["allele_type"].items():
                # increase allele type count
                sum_allele_type[type_of_allele] += 1
                # add allele name match to sample
                allele_match[allele] = (
                    type_of_allele + "_" + values["allele_match"][allele]
                )
            summary_result[sample] = sum_allele_type
            sample_allele_match[sample] = allele_match

    # save summary results to file
    with open(summary_result_file, "w") as fo:
        fo.write("Sample," + ",".join(allele_types) + "\n")
        for sample, counts in summary_result.items():
            fo.write(f"{sample},")
            for _, count in counts.items():
                fo.write(f"{count},")
            fo.write("\n")
    # save allele match to file
    with open(sample_allele_match_file, "w") as fo:
        fo.write("Sample," + ",".join(allele_list) + "\n")
        for sample, allele_cod in sample_allele_match.items():
            fo.write(f"{sample},")
            for allele in allele_list:
                fo.write(f"{allele_cod[allele]},")
            fo.write("\n")

    with open(sample_allele_detail_file, "w") as fo:
        fo.write(",".join(detail_heading) + "\n")
        for result in results:
            for sample, values in result.items():
                for allele, detail_value in values["allele_details"].items():
                    if type(detail_value[0]) is list:
                        for detail in detail_value:
                            fo.write(",".join(detail) + "\n")
                    else:
                        fo.write(",".join(detail_value) + "\n")
    if snp_request:
        snp_file = os.path.join(output, "snp_data.csv")
        with open(snp_file, "w") as fo:
            fo.write("Sample name,Locus name,Reference allele,Position,Base,Ref\n")
            for sample, values in result.items():
                # pdb.set_trace()
                for allele, snp_data in values["snp_data"].items():
                    for ref_allele, snp_info_list in snp_data.items():
                        # pdb.set_trace()
                        for snp_info in snp_info_list:
                            fo.write(
                                sample
                                + ","
                                + allele
                                + ","
                                + ref_allele
                                + ","
                                + ",".join(snp_info)
                                + "\n"
                            )
