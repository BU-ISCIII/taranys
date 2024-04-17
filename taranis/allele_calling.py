import io
import concurrent.futures
import logging
import os
import rich.console

import taranis.utils
import taranis.blast

from collections import OrderedDict
from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO

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
        hit_lenght_perc: float,
        perc_identity: int,
        out_folder: str,
        inf_alle_obj: object,
        snp_request: bool = False,
        aligment_request: bool = False,
        tpr_limit: int = 80,
        increase_sequence: int = 20,
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
            tpr_limit (int, optional): lower threshold to consider trunked proteine. Defaults to 80.
            increase_sequence (int, optional): increase sequence to be analysed. Defaults to 20.
        """
        self.prediction_data = annotation  # store prediction annotation
        self.sample_file = sample_file
        self.sample_contigs = taranis.utils.read_fasta_file(
            self.sample_file, convert_to_dict=True
        )
        self.schema = schema
        self.ref_alleles = reference_alleles
        self.hit_lenght_perc = hit_lenght_perc
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
        self.tpr_limit = tpr_limit / 100
        self.increase_sequence = increase_sequence

    def assign_allele_type(
        self,
        valid_blast_results: list,
        locus_file: str,
        locus_name: str,
        ref_allele_seq: str,
    ) -> list:
        """Assign allele type to the allele

        Args:
            valid_blast_results (list): information collected by running blast
            locus_file (str): file name with locus alleles sequences
            locus_name (str): locus name
            ref_allele_seq (str): reference allele sequence

        Returns:
            list: containing allele classification, allele match id, and allele
            details
        """

        def _check_plot(allele_details: dict) -> bool:
            """Check if allele is partial length

            Args:
                allele_details (dirt): allele details obtained with _get_allele_details() function. 

            Returns:
                bool: True if alignment is partial due to end of contig
            """
            if (
                allele_details["align_contig_start"] == "1"  # check  at contig start
                # check if contig ends is the same as match allele ends
                or allele_details["align_contig_end"] == allele_details["contig_length"]
                or allele_details["align_contig_end"] == "1"  # check reverse at contig end
                # check if contig start is the same as match allele start reverse
                or allele_details["align_contig_start"] == allele_details["contig_length"]
            ):
                return True
            return False

        def _extend_seq_find_start_stop_codon(
            direction: str,
            contig_seq: str,
            start: int,
            end: int,
            limit: int,
            search: str = "5_prime",
        ) -> list:
            """Extend match sequence, according to increase_sequence in order to try to
                find the stop or start codon.
            Args:

            Returns:

            """
            protein = "-"
            error = False
            error_details = "-"
            i = 0

            # Extend the sequence to find a valid start or stop codon
            if direction == "reverse":
                contig_seq = Seq(contig_seq).reverse_complement()
                start, end = len(contig_seq) - end, len(contig_seq) - start
            for _ in range(limit):
                i += 3
                if search == "5_prime":
                    extended_start = max(0, start - i)
                    extended_end = end
                elif search == "3_prime":
                    extended_start = start
                    extended_end = min(len(contig_seq), end + i)

                extended_seq = contig_seq[extended_start:extended_end]
                _, protein, error, error_details = taranis.utils.convert_to_protein(
                    extended_seq, force_coding=True
                )
                if not error:
                    return protein, extended_seq, extended_start, extended_end, error, error_details

            return protein, contig_seq[start:end], start, end, error, error_details

        def _get_allele_details(
            blast_result: str, locus_name: str, ref_allele_seq
        ) -> dict:
            """Collect blast details, add gene annotation, and protein sequence.

            Args:
                blast_result (str): information collected by running blast
                locus_name (str):  allele name
                ref_allele_seq (str): reference allele sequence

            Returns:
                dict:
                allele_details{
                        "sample_name": str,
                        "contig_name": str,
                        "locus_name": str,
                        "ref_allele_name": str,
                        "allele_type": str,
                        "ref_allele_length": str,
                        "alignment_length": str,
                        "contig_length": str,
                        "align_contig_start": str,
                        "align_contig_end": str,
                        "strand": str,
                        "sample_allele_seq": str,
                        "ref_allele_seq": str,
                        "gene_annotation": str,
                        "product_annot": str,
                        "ref_allele_quality": str,
                        "protein_seq": str,
                        "prot_strand": str,
                        "prot_error": bool,
                        "prot_error_details": str,
                    }
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

            if int(split_blast_result[10]) > int(split_blast_result[9]):
                strand = "+"
            else:
                strand = "-"

            # remove the gaps in sequences
            match_sequence = split_blast_result[13].replace("-", "")
            # check if the sequence is coding
            direction, protein, prot_error, prot_error_details = (
                taranis.utils.convert_to_protein(match_sequence, force_coding=True)
            )
            # get blast details
            allele_details = OrderedDict({
                "sample_name": self.s_name,  # sample name
                "locus_name": locus_name,  # core gene name
                "allele_type": "-",
                "ref_allele_name": split_blast_result[0],
                "contig_name": split_blast_result[1],  # contig name
                "contig_length": split_blast_result[15],
                "ref_allele_length": split_blast_result[3],
                "alignment_length": split_blast_result[4],
                "align_contig_start": split_blast_result[9],
                "align_contig_end": split_blast_result[10],
                "strand": strand,
                "sample_allele_seq": match_sequence,
                "ref_allele_seq": ref_allele_seq,
                "gene_annotation": gene_annotation,
                "product_annot": product_annotation,
                "ref_allele_quality": allele_quality,
                "protein_seq": protein,
                "prot_strand": direction,
                "prot_error": prot_error,
                "prot_error_details": prot_error_details,
            })
            return allele_details

        def _classify_allele(allele_file: str, match_sequence: str) -> str:
            """Find the allele name in the schema that match the sequence

            Args:
                allele_file (str): file with allele sequences
                match_sequence (str): sequence to be matched

            Returns:
                str: allele name in the schema that match the sequence
            """
            # Read the fasta file and create a dictionary mapping sequences to their record IDs
            sequence_dict = {
                str(record.seq): record.id
                for record in SeqIO.parse(allele_file, "fasta")
            }

            # Check if the match_sequence is in the dictionary and return the corresponding record ID part
            if match_sequence in sequence_dict:
                return sequence_dict[match_sequence]

            return False  # Return an empty string if no match is found

        match_allele_schema = ""

        if len(valid_blast_results) > 1:
            # could  be NIPHEM or NIPH
            sample_allele_data = []
            match_allele_seq = []

            for valid_blast_result in valid_blast_results:
                multi_allele_data = _get_allele_details(
                    valid_blast_result, locus_name, ref_allele_seq
                )
                # get match allele sequence
                match_allele_seq.append(multi_allele_data["sample_allele_seq"])
                sample_allele_data.append(multi_allele_data)
            if len(set(match_allele_seq)) == 1:
                # all sequuences are equal labelled as NIPHEM
                classification = "NIPHEM"
            else:
                # some of the sequences are different labelled as NIPH
                classification = "NIPH"
            # update coding allele type
            for idx in range(len(sample_allele_data)):
                sample_allele_data[idx]["allele_type"] = classification

        else:
            sample_allele_data = _get_allele_details(
                valid_blast_results[0], locus_name, ref_allele_seq
            )
            # found the allele in schema with the match sequence in the contig
            match_allele_schema = _classify_allele(
                locus_file, sample_allele_data["sample_allele_seq"]
            )
            # PLOT, TPR, ASM, ALM, INF, EXC are possible classifications
            if match_allele_schema:
                # exact match found labelled as EXC
                classification = "EXC"
                sample_allele_data["allele_type"] = classification + "_" + match_allele_schema
                return [
                    classification,
                    classification + "_" + match_allele_schema,
                    sample_allele_data,
                ]
            elif _check_plot(sample_allele_data):
                # match allele is partial length labelled as PLOT
                classification = "PLOT"
                sample_allele_data["allele_type"] = classification
                return [
                    classification,
                    classification,
                    sample_allele_data,
                ]

            if sample_allele_data["prot_error"]:
                # 0/1-based
                if sample_allele_data["strand"] == "+":
                    start = int(sample_allele_data["align_contig_start"]) - 1
                    end = int(sample_allele_data["align_contig_end"])
                else:
                    start = int(sample_allele_data["align_contig_end"]) - 1
                    end = int(sample_allele_data["align_contig_start"])

                # If strand "-" contig seq is reverse complemented but match sequence (split_blast_result[13])
                # is forward, so we need to change the contig seq to reverse_complement accordingly
                if sample_allele_data["strand"] == "-" and sample_allele_data["prot_strand"] == "reverse":
                    direction_contig = "forward"
                elif sample_allele_data["strand"] == "-" and sample_allele_data["prot_strand"] == "forward":
                    direction_contig = "reverse"
                else:
                    direction_contig = sample_allele_data["prot_strand"]

                if "is not a stop codon" in sample_allele_data["prot_error_details"]:
                    protein, extended_seq, new_start, new_end, prot_error, prot_error_details = (
                        _extend_seq_find_start_stop_codon(
                            direction=direction_contig,
                            contig_seq=self.sample_contigs[sample_allele_data["contig_name"]],
                            start=start,
                            end=end,
                            limit=self.increase_sequence,
                            search="3_prime",
                        )
                    )
                    if not prot_error:
                        sample_allele_data["protein_seq"] = protein
                        sample_allele_data["sample_allele_seq"] = extended_seq
                        sample_allele_data["align_contig_start"] = new_start
                        sample_allele_data["align_contig_end"] = new_end
                        sample_allele_data["prot_error"] = prot_error
                        sample_allele_data["prot_error_details"] = prot_error_details

                elif "is not a start codon" in sample_allele_data["prot_error_details"]:
                    protein, extended_seq, new_start, new_end, prot_error, prot_error_details = (
                        _extend_seq_find_start_stop_codon(
                            direction=direction_contig,
                            contig_seq=self.sample_contigs[sample_allele_data["contig_name"]],
                            start=start,
                            end=end,
                            limit=self.increase_sequence,
                            search="5_prime",
                        )
                    )
                    if not prot_error:
                        sample_allele_data["protein_seq"] = protein
                        sample_allele_data["sample_allele_seq"] = extended_seq
                        sample_allele_data["align_contig_start"] = new_start
                        sample_allele_data["align_contig_end"] = new_end
                        sample_allele_data["prot_error"] = prot_error
                        sample_allele_data["prot_error_details"] = prot_error_details

            if sample_allele_data["prot_error"]:
                classification = "TPR"
            # check if match allele is shorter than reference allele
            elif (
                int(len(sample_allele_data["sample_allele_seq"]))
                < int(sample_allele_data["ref_allele_length"]) - int(sample_allele_data["ref_allele_length"]) * 0.20
            ):
                classification = "ASM"
            # check if match allele is longer than reference allele
            elif (
                int(len(sample_allele_data["sample_allele_seq"]))
                > int(sample_allele_data["ref_allele_length"]) + int(sample_allele_data["ref_allele_length"]) * 0.20
            ):
                classification = "ALM"
            else:
                # if sequence was not found after running grep labelled as INF
                classification = "INF"
            # assign an identification value to the new allele
            if not match_allele_schema:
                match_allele_schema = str(
                    self.inf_alle_obj.get_inferred_allele(sample_allele_data["sample_allele_seq"], locus_name)
                )

            sample_allele_data["allele_type"] = classification + "_" + match_allele_schema

        return [
            classification,
            classification + "_" + match_allele_schema,
            sample_allele_data,
        ]

    def discard_low_threshold_results(self, blast_results: list) -> list:
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
            if (int(blast_split[4]) / int(blast_split[3])) >= self.hit_lenght_perc:
                valid_blast_result.append(b_result)
        return valid_blast_result

    def search_allele(self):
        """ Search reference allele in contig files and classify

        Args:
            self

        Returns:
            result = {
                "allele_type": {},
                "allele_match": {},
                "allele_details": {},
                "snp_data": {},
                "alignment_data": {},
            }
        """
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
                f"Processing allele {ref_allele}: {count} of {len(self.ref_alleles)}"
            )

            alleles = taranis.utils.read_fasta_file(ref_allele, convert_to_dict=True)
            match_found = False
            count_2 = 0
            for r_id, r_seq in alleles.items():
                count_2 += 1

                log.debug(f"Running blast for {count_2} of  {len(alleles)}")
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
                    valid_blast_results = self.discard_low_threshold_results(
                        blast_result
                    )
                    if len(valid_blast_results) > 0:
                        match_found = True
                        break
                # Close object and discard memory buffer
                query_file.close()

            locus_file = os.path.join(self.schema, os.path.basename(ref_allele))
            locus_name = Path(locus_file).stem

            if match_found:
                (
                    result["allele_type"][locus_name],
                    result["allele_match"][locus_name],
                    result["allele_details"][locus_name],
                ) = self.assign_allele_type(
                    valid_blast_results, locus_file, locus_name, r_seq
                )
            else:
                # Sample does not have a reference allele to be matched
                # Keep LNF info
                result["allele_type"][locus_name] = "LNF"
                result["allele_match"][locus_name] = "LNF"
                details = OrderedDict()
                details["sample_name"] = self.s_name
                details["locus_name"] = locus_name
                details["allele_type"] = "LNF"
                result["allele_details"][locus_name] = details

            # prepare the data for snp and alignment analysis
            if result["allele_type"][locus_name] not in [
                "PLOT",
                "LNF",
                "NIPH",
                "NIPHEM",
            ]:
                try:
                    ref_allele_seq = result["allele_details"][locus_name]["ref_allele_seq"]
                except KeyError as e:
                    log.error("Error in allele details")
                    log.error(e)
                    stderr.print(f"Error in allele details{e}")
                    continue

                allele_seq = result["allele_details"][locus_name]["sample_allele_seq"]
                ref_allele_name = result["allele_details"][locus_name]["ref_allele_name"]

                if self.snp_request and result["allele_type"][locus_name] != "LNF":
                    # run snp analysis
                    result["snp_data"][locus_name] = taranis.utils.get_snp_information(
                        ref_allele_seq, allele_seq, ref_allele_name
                    )

                if self.aligment_request and result["allele_type"][locus_name] != "LNF":
                    # run alignment analysis
                    result["alignment_data"][locus_name] = (
                        taranis.utils.get_alignment_data(
                            ref_allele_seq, allele_seq, ref_allele_name
                        )
                    )
        # delete blast folder
        _ = taranis.utils.delete_folder(os.path.join(self.blast_dir, self.s_name))
        return result


def parallel_execution(
    sample_file: str,
    schema: str,
    prediction_data: dict,
    reference_alleles: list,
    hit_lenght_perc: float,
    perc_identity: int,
    out_folder: str,
    inf_alle_obj: object,
    snp_request: bool = False,
    aligment_request: bool = False,
    trp_limit: int = 80,
    increase_sequence: int = 20,
):
    allele_obj = AlleleCalling(
        sample_file,
        schema,
        prediction_data,
        reference_alleles,
        hit_lenght_perc,
        perc_identity,
        out_folder,
        inf_alle_obj,
        snp_request,
        aligment_request,
        trp_limit,
        increase_sequence,
    )
    sample_name = Path(sample_file).stem
    stderr.print(f"[green] Analyzing sample {sample_name}")
    log.info(f"Analyzing sample {sample_name}")
    return {sample_name: allele_obj.search_allele()}


def create_multiple_alignment(
    ref_alleles_seq: dict, results: list, a_list: str, alignment_folder: str, mafft_cpus
) -> None:
    allele_multiple_align = []
    for ref_id, ref_seq in ref_alleles_seq[a_list].items():
        input_buffer = StringIO()
        # get the reference allele sequence
        input_buffer.write(">Ref_" + ref_id + "\n")
        input_buffer.write(ref_seq + "\n")
        # get the sequences for sample on the same allele
        for result in results:
            for sample, values in result.items():
                # discard the allele if it is LNF
                if values["allele_type"][a_list] == "LNF":
                    continue
                # get the allele name in sample
                input_buffer.write(
                    ">"
                    + sample
                    + "_"
                    + a_list
                    + "_"
                    + values["allele_details"][a_list][4]
                    + "\n"
                )
                # get the sequence of the allele in sample
                input_buffer.write(values["allele_details"][a_list][15] + "\n")
        # print(input_buffer.tell())
        input_buffer.seek(0)

        allele_multiple_align.append(
            taranis.utils.get_multiple_alignment(input_buffer, mafft_cpus)
        )
        # release memory
        input_buffer.close()
    # save multiple alignment to file
    with open(
        os.path.join(alignment_folder, a_list + "_multiple_alignment.aln"), "w"
    ) as fo:
        for alignment in allele_multiple_align:
            for align in alignment:
                fo.write(align)


def collect_data(
    results: list,
    output: str,
    snp_request: bool,
    aligment_request: bool,
    ref_alleles: list,
    cpus: int,
) -> None:
    """Collect data for the allele calling analysis, done for each sample and
    create the summary file, graphics, and if requested snp and alignment files

    Args:
        results (list): list of allele calling data results for each sample
        output (str): output folder
        snp_request (bool): request to save snp to file
        aligment_request (bool): request to save alignment and multi alignemte to file
        ref_alleles (list): reference alleles
        cpus (int): number of cpus to be used if alignment is requested
    """

    def stats_graphics(stats_folder: str, summary_result: dict) -> None:
        stderr.print("Creating graphics")
        log.info("Creating graphics")
        allele_types = [
            "NIPHEM",
            "NIPH",
            "EXC",
            "PLOT",
            "ASM",
            "ALM",
            "INF",
            "LNF",
            "TPR",
        ]
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

    def read_reference_alleles(ref_alleles: list) -> dict[dict]:
        # read reference alleles
        ref_alleles_data = {}
        for ref_allele in ref_alleles:
            alleles = {}
            with open(ref_allele, "r") as fh:
                for record in SeqIO.parse(fh, "fasta"):
                    alleles[record.id] = str(record.seq)
            ref_alleles_data[Path(ref_allele).stem] = alleles
        return ref_alleles_data

    summary_result_file = os.path.join(output, "allele_calling_summary.csv")
    allele_matrix_file = os.path.join(output, "allele_calling_match.csv")
    allele_detail_file = os.path.join(output, "contig_alignment_info.csv")
    allele_types = ["NIPHEM", "NIPH", "EXC", "PLOT", "ASM", "ALM", "INF", "LNF", "TPR"]
    allele_detail_heading = [
        "sample_name",
        "locus_name",
        "allele_type",
        "ref_allele_name",
        "contig_name",
        "contig_length",
        "ref_allele_length",
        "alignment_length",
        "align_contig_start",
        "align_contig_end",
        "strand",
        "sample_allele_seq",
        "ref_allele_seq",
        "gene_annotation",
        "product_annot",
        "ref_allele_quality",
        "protein_seq",
        "prot_strand",
        "prot_error",
        "prot_error_details",
    ]

    summary_result = {}  # used for summary file and allele classification graphics
    allele_matrix_result = {}  # used for allele match file

    # get allele list
    locus_list = [Path(ref_allele).stem for ref_allele in ref_alleles]

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
                    values["allele_match"][allele]
                )
            summary_result[sample] = sum_allele_type
            allele_matrix_result[sample] = allele_match

    # save summary results to file
    with open(summary_result_file, "w") as fo:
        fo.write("Sample," + ",".join(allele_types) + "\n")
        for sample, counts in summary_result.items():
            fo.write(f"{sample},")
            for _, count in counts.items():
                fo.write(f"{count},")
            fo.write("\n")

    # save allele match to file
    with open(allele_matrix_file, "w") as fo:
        fo.write("Sample," + ",".join(locus_list) + "\n")
        for sample, allele_cod in allele_matrix_result.items():
            fo.write(f"{sample}")
            for allele in locus_list:
                fo.write(f",{allele_cod[allele]}")
            fo.write("\n")

    with open(allele_detail_file, "w") as fo:
        fo.write(",".join(allele_detail_heading) + "\n")
        for result in results:
            for sample, values in result.items():
                for allele, detail_value in values["allele_details"].items():
                    if type(detail_value) is list:
                        for detail in detail_value:
                            if detail["allele_type"] != "LNF":
                                fo.write(",".join([str(value) for value in detail.values()]) + "\n")
                    else:
                        if detail_value["allele_type"] != "LNF":
                            fo.write(",".join([str(value) for value in detail_value.values()]) + "\n")

    # save snp to file if requested
    if snp_request:
        for result in results:
            for sample, values in result.items():
                snp_file = os.path.join(output, sample + "_snp_data.csv")
                with open(snp_file, "w") as fo:
                    fo.write(
                        "Sample name,Locus name,Reference allele,Position,Ref,Alt,Codon Ref,Codon Alt,Amino Ref,Amino Alt,Category Ref,Category Alt\n"
                    )
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
        align_collection = {}
        for result in results:
            for sample, values in result.items():
                for allele, alignment_data in values["alignment_data"].items():
                    if allele not in align_collection:
                        align_collection[allele] = OrderedDict()

                    # align_collection[allele][sample] = []
                    for _, value in alignment_data.items():
                        align_collection[allele][sample] = value
        # save alignment to file
        for allele, samples in align_collection.items():
            with open(os.path.join(alignment_folder, allele + ".txt"), "w") as fo:
                for sample, alignment_data in samples.items():
                    fo.write(allele + "_sample_" + sample + "\n")
                    fo.write("\n".join(alignment_data) + "\n")

        # create multiple alignment files
        stderr.print("Processing multiple alignment information")
        log.info("Processing multiple alignment information")
        ref_alleles_seq = read_reference_alleles(ref_alleles)
        # assign cpus to be used in multiple alignment
        mul_align_cpus = 1 if cpus // 3 == 0 else cpus // 3
        mafft_cpus = 1 if mul_align_cpus == 1 else 3
        m_align = []
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=mul_align_cpus
        ) as executor:
            futures = [
                executor.submit(
                    create_multiple_alignment,
                    ref_alleles_seq,
                    results,
                    a_list,
                    alignment_folder,
                    mafft_cpus,
                )
                for a_list in locus_list
            ]
        for future in concurrent.futures.as_completed(futures):
            try:
                m_align.append(future.result())
            except Exception as e:
                print(e)
                continue

        # for a_list in allele_list:
        #     _ = create_multiple_alignment(ref_alleles_seq, results, a_list, alignment_folder)

    # Create graphics
    stats_graphics(output, summary_result)
    return
