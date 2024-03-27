import io
import logging
import os
import rich.console

import taranis.utils
import taranis.blast

from collections import OrderedDict
from pathlib import Path
from Bio import SeqIO
import pdb

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
        threshold: float,
        perc_identity: int,
        out_folder: str,
        inf_alle_obj: object,
        snp_request: bool = False,
        aligment_request: bool = False,
    ):
        """Allele calling initial creation object

        Args:
            sample_file (str): assembly file
            schema (str): folder with alleles schema
            annotation (dict): annotation of locus according to prokka
            reference_alleles (list): folder with reference alleles
            threshold (float): threshold to consider a match in blast
            out_folder (str): output folder
            inf_alle_obj (object): object to infer alleles
            snp_request (bool, optional): snp saved to file. Defaults to False.
            aligment_request (bool, optional): allignment saved to file. Defaults to False.
        """
        self.prediction_data = annotation  # store prediction annotation
        self.sample_file = sample_file
        self.schema = schema
        self.ref_alleles = reference_alleles
        self.threshold = threshold
        self.perc_identity = perc_identity
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
        self, blast_results: list, allele_file: str, allele_name: str
    ) -> list:
        """Assign allele type to the allele

        Args:
            blast_result (list): information collected by running blast
            allele_file (str): file name with allele sequence
            allele_name (str): allele name

        Returns:
            list: containing allele classification, allele match id, and allele
            details
        """

        def check_if_plot(column_blast_res: list) -> bool:
            """Check if allele is partial length

            Args:
                column_blast_res (list): blast result

            Returns:
                bool: True if allele is partial length
            """
            if (
                column_blast_res[9] == "1"  # check  at contig start
                # check if contig ends is the same as match allele ends
                or column_blast_res[15] == column_blast_res[10]
                or column_blast_res[10] == "1"  # check reverse at contig end
                # check if contig start is the same as match allele start reverse
                or column_blast_res[9] == column_blast_res[15]
            ):
                return True
            return False

        def discard_low_threshold_results(blast_results: list) -> list:
            """Discard blast results with lower threshold

            Args:
                blast_results (list): blast results

            Returns:
                list: blast results with higher query size
            """
            valid_blast_result = []
            for b_result in blast_results:
                blast_split = b_result.split("\t")
                # check if the division of the match contig length by the
                # reference allele length is higher than the threshold
                # pdb.set_trace()
                if (int(blast_split[4]) / int(blast_split[3])) >= self.threshold:
                    valid_blast_result.append(b_result)
            return valid_blast_result

        def get_blast_details(blast_result: str, allele_name: str) -> list:
            """Collect blast details and modify the order of the columns

            Args:
                blast_result (str): information collected by running blast
                allele_name (str):  allele name

            Returns:
                list: containing allele details in the correct order to be saved
                    blast_details[0] = sample name
                    blast_details[1] = contig name
                    blast_details[2] = core gene name
                    blast_details[3] = allele gene
                    blast_details[4] = coding allele type
                    blast_details[5] = reference allele length
                    blast_details[6] = match alignment length
                    blast_details[7] = contig length
                    blast_details[8] = match contig position start
                    blast_details[9] = match contig position end
                    blast_details[10] = direction
                    blast_details[11] = gene annotation
                    blast_details[12] = product annotation
                    blast_details[13] = allele quality
                    blast_details[14] = match sequence in contig
                    blast_details[15] = reference allele sequence
            """
            split_blast_result = blast_result.split("\t")
            match_allele_name = split_blast_result[0]
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
            # pdb.set_trace()
            if int(split_blast_result[10]) > int(split_blast_result[9]):
                direction = "+"
            else:
                direction = "-"
            # get blast details
            blast_details = [
                self.s_name,  # sample name
                split_blast_result[1],  # contig name
                allele_name,  # core gene name
                split_blast_result[0],  # allele gene
                "coding",  # coding allele type. To be filled later idx = 4
                split_blast_result[3],  # reference allele length
                split_blast_result[4],  # match alignment length
                split_blast_result[15],  # contig length
                split_blast_result[9],  # match contig position start
                split_blast_result[10],  # match contig position end
                direction,
                gene_annotation,
                product_annotation,
                allele_quality,
                split_blast_result[13],  # match sequence in contig
                split_blast_result[15],  # reference allele sequence
            ]
            # pdb.set_trace()
            return blast_details

        def find_match_allele_schema(allele_file: str, match_sequence: str) -> str:
            """Find the allele name in the schema that match the sequence

            Args:
                allele_file (str): file with allele sequences
                match_sequence (str): sequence to be matched

            Returns:
                str: allele name in the schema that match the sequence
            """
            grep_result = taranis.utils.grep_execution(
                allele_file, match_sequence, "-b1"
            )
            if len(grep_result) > 0:
                return grep_result[0].split("_")[1]
            return ""

        valid_blast_results = discard_low_threshold_results(blast_results)
        match_allele_schema = ""
        if len(valid_blast_results) == 0:
            # no match results labelled as LNF. details data filled with empty data
            return ["LNF", "-", ["-," * 15]]
        if len(valid_blast_results) > 1:
            # could  be NIPHEM or NIPH
            b_split_data = []
            match_allele_seq = []
            for valid_blast_result in valid_blast_results:
                multi_allele_data = get_blast_details(valid_blast_result, allele_name)
                # get match allele sequence
                match_allele_seq.append(multi_allele_data[14])
                b_split_data.append(multi_allele_data)
                # check if match allele is in schema
                if match_allele_schema == "":
                    # find the allele in schema with the match sequence in the contig
                    match_allele_schema = find_match_allele_schema(
                        allele_file, multi_allele_data[14]
                    )
            if len(set(match_allele_seq)) == 1:
                # all sequuences are equal labelled as NIPHEM
                classification = "NIPHEM"
            else:
                # some of the sequences are different labelled as NIPH
                classification = "NIPH"
            # update coding allele type
            for (idx,) in range(len(b_split_data)):
                b_split_data[idx][4] = classification + "_" + match_allele_schema
        else:
            b_split_data = get_blast_details(valid_blast_results[0], allele_name)
            # found the allele in schema with the match sequence in the contig
            match_allele_schema = find_match_allele_schema(
                allele_file, b_split_data[14]
            )
            # PLOT, ASM, ALM, INF, EXC are possible classifications
            if check_if_plot(b_split_data):
                # match allele is partial length labelled as PLOT
                classification = "PLOT"

            # check if match allele is shorter than reference allele
            elif int(b_split_data[5]) < int(b_split_data[6]):
                classification = "ASM"
            # check if match allele is longer than reference allele
            elif int(b_split_data[5]) > int(b_split_data[6]):
                classification = "ALM"
            else:
                # if sequence was not found after running grep labelled as INF
                if match_allele_schema == "":
                    classification = "INF"
                else:
                    # exact match found labelled as EXC
                    classification = "EXC"
            # assign an identification value to the new allele
            if match_allele_schema == "":
                match_allele_schema = str(
                    self.inf_alle_obj.get_inferred_allele(b_split_data[14], allele_name)
                )
        # pdb.set_trace()
        b_split_data[4] = classification + "_" + match_allele_schema
        return [
            classification,
            classification + "_" + match_allele_schema,
            b_split_data,
        ]

    def search_match_allele(self):
        # Create  blast db with sample file

        result = {
            "allele_type": {},
            "allele_match": {},
            "allele_details": {},
            "snp_data": {},
            "alignment_data": {},
        }
        count = 0
        for ref_allele in self.ref_alleles:
            count += 1
            log.debug(
                " Processing allele ",
                ref_allele,
                " ",
                count,
                " of ",
                len(self.ref_alleles),
            )

            alleles = OrderedDict()
            match_found = False
            with open(ref_allele, "r") as fh:
                for record in SeqIO.parse(fh, "fasta"):
                    alleles[record.id] = str(record.seq)
            count_2 = 0
            for r_id, r_seq in alleles.items():
                count_2 += 1

                log.debug("Running blast for ", count_2, " of ", len(alleles))
                # create file in memory to increase speed
                query_file = io.StringIO()
                query_file.write(">" + r_id + "\n" + r_seq)
                query_file.seek(0)
                blast_result = self.blast_obj.run_blast(
                    query_file.read(),
                    perc_identity=self.perc_identity,
                    num_threads=1,
                    query_type="stdin",
                )
                if len(blast_result) > 0:
                    match_found = True
                    break
            # Close object and discard memory buffer
            query_file.close()
            if match_found:
                allele_file = os.path.join(self.schema, os.path.basename(ref_allele))
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
            if self.aligment_request and result["allele_type"][allele_name] == "INF":
                # run alignment analysis
                allele_seq = result["allele_details"][allele_name][14]
                result["alignment_data"][
                    allele_name
                ] = taranis.utils.get_alignment_data(allele_seq, alleles)
        return result


def parallel_execution(
    sample_file: str,
    schema: str,
    prediction_data: dict,
    reference_alleles: list,
    threshold: float,
    perc_identity: int,
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
        threshold,
        perc_identity,
        out_folder,
        inf_alle_obj,
        snp_request,
        aligment_request,
    )
    sample_name = Path(sample_file).stem
    stderr.print(f"[green] Analyzing sample {sample_name}")
    log.info(f"Analyzing sample {sample_name}")
    return {sample_name: allele_obj.search_match_allele()}


def collect_data(
    results: list, output: str, snp_request: bool, aligment_request: bool
) -> None:
    def stats_graphics(stats_folder: str, summary_result: dict) -> None:
        stderr.print("Creating graphics")
        log.info("Creating graphics")
        allele_types = ["NIPHEM", "NIPH", "EXC", "PLOT", "ASM", "ALM", "INF", "LNF"]
        # inizialize classification data
        classif_data = {}
        for allele_type in allele_types:
            classif_data[allele_type] = []
        graphic_folder = os.path.join(stats_folder, "graphics")

        _ = taranis.utils.create_new_folder(graphic_folder)
        s_list = []
        # collecting data to create graphics
        for sample, classif_counts in summary_result.items():
            s_list.append(sample)  # create list of samples
            for classif, count in classif_counts.items():
                classif_data[classif].append(int(count))
        # create graphics per each classification type
        for allele_type, counts in classif_data.items():
            _ = taranis.utils.create_graphic(
                graphic_folder,
                str(allele_type + "_graphic.png"),
                "bar",
                s_list,
                counts,
                ["Samples", "number"],
                str("Number of " + allele_type + " in samples"),
            )
        return

    summary_result_file = os.path.join(output, "allele_calling_summary.csv")
    sample_allele_match_file = os.path.join(output, "allele_calling_match.csv")
    sample_allele_detail_file = os.path.join(output, "matching_contig.csv")
    allele_types = ["NIPHEM", "NIPH", "EXC", "PLOT", "ASM", "ALM", "INF", "LNF"]
    detail_heading = [
        "sample",
        "contig",
        "core gene",
        "reference allele name",
        "codification",
        "query length",
        "match length",
        "contig length",
        "contig start",
        "contig stop",
        "direction",
        "gene notation",
        "product notation",
        "allele quality",
        "sequence",
    ]

    summary_result = {}  # used for summary file and allele classification graphics
    sample_allele_match = {}  # used for allele match file

    # get allele list
    # pdb.set_trace()
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
                    # type_of_allele + "_" + values["allele_match"][allele]
                    values["allele_match"][allele]
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
            fo.write(f"{sample}")
            for allele in allele_list:
                fo.write(f",{allele_cod[allele]}")
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
                for allele, snp_data in values["snp_data"].items():
                    for ref_allele, snp_info_list in snp_data.items():
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
    # create alignment files
    if aligment_request:
        alignment_folder = os.path.join(output, "alignments")
        _ = taranis.utils.create_new_folder(alignment_folder)
        for result in results:
            for sample, values in result.items():
                for allele, alignment_data in values["alignment_data"].items():
                    with open(
                        os.path.join(alignment_folder, sample + "_" + allele + ".txt"),
                        "w",
                    ) as fo:
                        for ref_allele, alignments in alignment_data.items():
                            fo.write(ref_allele + "\n")
                            for alignment in alignments:
                                fo.write(alignment + "\n")

    # Create graphics
    stats_graphics(output, summary_result)
