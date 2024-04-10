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
        threshold: float,
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
        self.sample_records = taranis.utils.read_fasta_file(
            self.sample_file, convert_to_dict=True
        )
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
        self.tpr_limit = tpr_limit / 100
        self.increase_sequence = increase_sequence * 3

    def assign_allele_type(
        self,
        valid_blast_results: list,
        allele_file: str,
        allele_name: str,
        ref_allele_seq: str,
    ) -> list:
        """Assign allele type to the allele

        Args:
            valid_blast_results (list): information collected by running blast
            allele_file (str): file name with allele sequence
            allele_name (str): allele name
            ref_allele_seq (str): reference allele sequence

        Returns:
            list: containing allele classification, allele match id, and allele
            details
        """

        def _check_if_plot(column_blast_res: list) -> bool:
            """Check if allele is partial length

            Args:
                column_blast_res (list): blast result

            Returns:
                bool: True if allele is partial length
            """
            if (
                column_blast_res[8] == "1"  # check  at contig start
                # check if contig ends is the same as match allele ends
                or column_blast_res[9] == column_blast_res[7]
                or column_blast_res[9] == "1"  # check reverse at contig end
                # check if contig start is the same as match allele start reverse
                or column_blast_res[8] == column_blast_res[7]
            ):
                return True
            return False

        def _extend_sequence_for_finding_start_stop_codon(
            split_blast_result: list,
            prot_error_result: str,
            predicted_prot_seq: str,
            search_codon: str = "stop",
        ) -> list:
            """Extend match sequence, according the (increase_sequence) for
                trying find the stop or start codon. When parameter is set to
                stop additional nucleotides are added to extend the chance to
                find out the codon stop.
                If parameter is set to start then additional nucleotide is added
                on the start value to identify that is a valid start codon. If
                true then additional nucletotides are added to find the stop codon.

            Args:
                split_blast_result (list): list having the informaction collected
                    by running blast
                prot_error_result (str): protein conversion result
                predicted_prot_seq (str): predicted protein sequence
                search_codon (str, optional): codon to be found. 2 values are
                    allowed start of stop. By default is stop.

            Returns:
                list: updated information if stop or start codon is found and the
                updated protein sequence and protein conversion result if changed
            """
            # collect data for checking PLOT
            data_for_plot = [""] * 10
            # cop the contig length
            data_for_plot[7] = split_blast_result[15]
            # copy start position
            data_for_plot[8] = split_blast_result[9]
            # copy end position
            data_for_plot[9] = split_blast_result[10]
            # check if PLOT
            if not _check_if_plot(data_for_plot):
                # remove the "-" character in the contig sequence in case that
                # there was not possible to find the start/stop codon and
                # function return the original blast result
                split_blast_result[13] = split_blast_result[13].replace("-", "")
                # fetch the sequence until the last triplet is stop codon
                contig_seq = self.sample_records[split_blast_result[1]]
                start_seq = int(split_blast_result[9])
                stop_seq = int(split_blast_result[10])
                if stop_seq > start_seq:
                    # sequence direction is forward
                    direction = "forward"
                    if search_codon == "start":
                        # add nucleotides according to the first match in the
                        # reference allele
                        start_ref_allele = int(split_blast_result[11]) // 3 * 3
                        # try extended 1 nucleotide to find the start codon
                        start_seq_found = False
                        # subtract 1 because index start at 0
                        new_start_seq = start_seq - start_ref_allele - 1
                        for _ in range(1, 3):
                            new_start_seq -= 1
                            if (
                                contig_seq[new_start_seq : new_start_seq + 3]
                                in taranis.utils.START_CODON_FORWARD
                            ):
                                # increase 1 because we substact 1 when searching
                                # for stop codon
                                start_seq = new_start_seq + 1
                                start_seq_found = True
                                break
                                # continue to find the stop codon with the new start
                        if not start_seq_found:
                            # start codon not found. Return the original blast result
                            return (
                                split_blast_result,
                                prot_error_result,
                                predicted_prot_seq,
                            )

                    # adjust the sequence to be a triplet
                    interval = (stop_seq - start_seq) // 3 * 3
                    new_stop_seq = start_seq + interval + self.increase_sequence
                    start_seq -= 1
                    # if the increased length is higher than the contig length
                    # adjust the stop sequence to maximun contig length
                    # multiply by 3.
                    if stop_seq > len(contig_seq):
                        stop_seq = len(contig_seq) // 3 * 3
                    else:
                        stop_seq = new_stop_seq - 1
                    c_sequence = contig_seq[start_seq:stop_seq]
                else:
                    # sequence direction is reverse
                    direction = "reverse"
                    if search_codon == "start":
                        if (
                            contig_seq[start_seq - 2 : start_seq + 1]
                            in taranis.utils.START_CODON_REVERSE
                        ):
                            start_seq += 1
                            # continue to find the stop codon with the new start
                        else:
                            # start codon not found. Return the original blast result
                            return (
                                split_blast_result,
                                prot_error_result,
                                predicted_prot_seq,
                            )
                    # adjust the sequence to be a triplet
                    interval = (start_seq - stop_seq) // 3 * 3
                    new_stop_seq = start_seq - interval - self.increase_sequence
                    # if the increased length is lower than 0 (contig start)
                    # position, adjust the start sequence to minumum contig
                    # length multiply by 3
                    if new_stop_seq < 0:
                        # get the minimum contig length that is multiple by 3
                        stop_seq = stop_seq % 3 - 1
                    else:
                        stop_seq = new_stop_seq
                    # get the sequence in reverse
                    c_sequence = str(
                        Seq(contig_seq[stop_seq:start_seq]).reverse_complement()
                    )
                new_prot_conv_result = taranis.utils.convert_to_protein(
                    c_sequence, force_coding=False, delete_incompleted_triplet=False
                )
                # check if stop codon is found in protein sequence

                if (
                    "protein" in new_prot_conv_result
                    and "*" in new_prot_conv_result["protein"]
                ):
                    # increase 3 nucleotides beecause index start at 0
                    new_seq_length = new_prot_conv_result["protein"].index("*") * 3 + 3
                    match_sequence = c_sequence[:new_seq_length]
                    split_blast_result[4] = str(new_seq_length)
                    split_blast_result[13] = match_sequence
                    prot_error_result = "-"
                    predicted_prot_seq = new_prot_conv_result["protein"][
                        0 : new_seq_length // 3
                    ]
                    # update the start and stop position
                    if direction == "forward":
                        split_blast_result[10] = str(
                            int(split_blast_result[9]) + new_seq_length
                        )
                    else:
                        split_blast_result[10] = str(
                            int(split_blast_result[9]) - new_seq_length
                        )
                # ignore the previous process if stop codon is not found

            return split_blast_result, prot_error_result, predicted_prot_seq

        def _get_blast_details(
            blast_result: str, allele_name: str, ref_allele_seq
        ) -> list:
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
                    blast_details[14] = protein conversion result
                    blast_details[15] = match sequence in contig
                    blast_details[16] = reference allele sequence
                    blast_details[17] = predicted protein sequence
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
                direction = "+"
            else:
                direction = "-"
            # remove the gaps in sequences
            match_sequence = split_blast_result[13].replace("-", "")
            # check if the sequence is coding
            prot_conv_result = taranis.utils.convert_to_protein(
                match_sequence, force_coding=False, delete_incompleted_triplet=True
            )
            prot_error_result = (
                prot_conv_result["error"] if "error" in prot_conv_result else "-"
            )
            predicted_prot_seq = (
                prot_conv_result["protein"] if "protein" in prot_conv_result else "-"
            )
            # remove if extra nucleotides are added at the end of the last
            # completed triplet before the stop codon
            if "extra nucleotides after stop codon" in prot_error_result:
                new_seq_len = len(match_sequence) // 3 * 3
                match_sequence = match_sequence[:new_seq_len]
                split_blast_result[4] = str(new_seq_len)
                # reset the error message
                prot_error_result = "-"

            # extend the sequence to find the stop codon
            elif "Last triplet sequence is not a stop codon" in prot_error_result:
                (
                    split_blast_result,
                    prot_error_result,
                    predicted_prot_seq,
                ) = _extend_sequence_for_finding_start_stop_codon(
                    split_blast_result,
                    prot_error_result,
                    predicted_prot_seq,
                    search_codon="stop",
                )
                # update the match sequence
                match_sequence = split_blast_result[13]
            # extend the sequence to find the start codon
            elif "Sequence does not have a start codon" in prot_error_result:
                (
                    split_blast_result,
                    prot_error_result,
                    predicted_prot_seq,
                ) = _extend_sequence_for_finding_start_stop_codon(
                    split_blast_result,
                    prot_error_result,
                    predicted_prot_seq,
                    search_codon="start",
                )
                # update the match sequence
                match_sequence = split_blast_result[13]
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
                prot_error_result,  # protein conversion result
                match_sequence,  # match sequence in contig
                ref_allele_seq,  # reference allele sequence
                predicted_prot_seq,  # predicted protein sequence
            ]
            return blast_details

        def _find_match_allele_schema(allele_file: str, match_sequence: str) -> str:
            """Find the allele name in the schema that match the sequence

            Args:
                allele_file (str): file with allele sequences
                match_sequence (str): sequence to be matched

            Returns:
                str: allele name in the schema that match the sequence
            """
            grep_result = taranis.utils.grep_execution(
                allele_file, match_sequence, "-xb1"
            )
            if len(grep_result) > 0:
                return grep_result[0].split("_")[1]
            return ""

        # valid_blast_results = _discard_low_threshold_results(blast_results)
        match_allele_schema = ""
        # if len(valid_blast_results) == 0:
        # no match results labelled as LNF. details data filled with empty data
        #    return ["LNF", "LNF", ["-"] * 18]
        if len(valid_blast_results) > 1:
            # could  be NIPHEM or NIPH
            b_split_data = []
            match_allele_seq = []
            for valid_blast_result in valid_blast_results:
                multi_allele_data = _get_blast_details(
                    valid_blast_result, allele_name, ref_allele_seq
                )
                # get match allele sequence
                match_allele_seq.append(multi_allele_data[14])
                b_split_data.append(multi_allele_data)
                # check if match allele is in schema
                if match_allele_schema == "":
                    # find the allele in schema with the match sequence in the contig
                    match_allele_schema = _find_match_allele_schema(
                        allele_file, multi_allele_data[15]
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
            b_split_data = _get_blast_details(
                valid_blast_results[0], allele_name, ref_allele_seq
            )
            # found the allele in schema with the match sequence in the contig
            match_allele_schema = _find_match_allele_schema(
                allele_file, b_split_data[15]
            )

            # PLOT, TPR, ASM, ALM, INF, EXC are possible classifications
            if match_allele_schema != "":
                # exact match found labelled as EXC
                classification = "EXC"
            elif _check_if_plot(b_split_data):
                # match allele is partial length labelled as PLOT
                classification = "PLOT"
            # check if protein length divided by the length of triplet matched
            # sequence is lower the the tpr limit
            elif (
                b_split_data[14] == "Multiple stop codons"
                and b_split_data[17].index("*") / (int(b_split_data[6]) / 3)
                < self.tpr_limit
            ):
                # labelled as TPR
                classification = "TPR"
                # check if match allele is shorter than reference allele
            elif int(b_split_data[6]) < int(b_split_data[5]):
                classification = "ASM"
            # check if match allele is longer than reference allele
            elif (
                int(b_split_data[6]) > int(b_split_data[5])
                or b_split_data[14] == "Last sequence is not a stop codon"
            ):
                classification = "ALM"
            else:
                # if sequence was not found after running grep labelled as INF
                classification = "INF"

            # assign an identification value to the new allele
            if match_allele_schema == "":
                match_allele_schema = str(
                    self.inf_alle_obj.get_inferred_allele(b_split_data[14], allele_name)
                )
        b_split_data[4] = classification + "_" + match_allele_schema
        return [
            classification,
            classification + "_" + match_allele_schema,
            b_split_data,
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
            if (int(blast_split[4]) / int(blast_split[3])) >= self.threshold:
                valid_blast_result.append(b_result)
        return valid_blast_result

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

            alleles = taranis.utils.read_fasta_file(ref_allele, convert_to_dict=True)
            match_found = False
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
                    valid_blast_results = self.discard_low_threshold_results(
                        blast_result
                    )
                    if len(valid_blast_results) > 0:
                        match_found = True
                        break
            # Close object and discard memory buffer
            query_file.close()
            allele_file = os.path.join(self.schema, os.path.basename(ref_allele))
            allele_name = Path(allele_file).stem
            if match_found:
                (
                    result["allele_type"][allele_name],
                    result["allele_match"][allele_name],
                    result["allele_details"][allele_name],
                ) = self.assign_allele_type(
                    valid_blast_results, allele_file, allele_name, r_seq
                )
            else:
                # Sample does not have a reference allele to be matched
                # Keep LNF info
                result["allele_type"][allele_name] = "LNF"
                result["allele_match"][allele_name] = allele_name
                details = ["-"] * 18
                details[0] = self.s_name
                details[2] = allele_name
                details[4] = "LNF"
                result["allele_details"][allele_name] = details
            # prepare the data for snp and alignment analysis
            try:
                ref_allele_seq = result["allele_details"][allele_name][16]
            except KeyError as e:
                log.error("Error in allele details")
                log.error(e)
                stderr.print(f"Error in allele details{e}")
                continue
            allele_seq = result["allele_details"][allele_name][15]
            ref_allele_name = result["allele_details"][allele_name][3]

            if self.snp_request and result["allele_type"][allele_name] != "LNF":
                # run snp analysis
                result["snp_data"][allele_name] = taranis.utils.get_snp_information(
                    ref_allele_seq, allele_seq, ref_allele_name
                )
            if self.aligment_request and result["allele_type"][allele_name] != "LNF":
                # run alignment analysis
                result["alignment_data"][allele_name] = (
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
    threshold: float,
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
        threshold,
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
    return {sample_name: allele_obj.search_match_allele()}


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
    sample_allele_match_file = os.path.join(output, "allele_calling_match.csv")
    sample_allele_detail_file = os.path.join(output, "matching_contig.csv")
    allele_types = ["NIPHEM", "NIPH", "EXC", "PLOT", "ASM", "ALM", "INF", "LNF", "TPR"]
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
        "reference allele quality",
        "protein conversion result",
        "match sequence",
        "reference allele sequence",
        "predicted protein sequence",
    ]

    summary_result = {}  # used for summary file and allele classification graphics
    sample_allele_match = {}  # used for allele match file

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
                for a_list in allele_list
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
