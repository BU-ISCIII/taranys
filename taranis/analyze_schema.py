import logging
import pandas as pd
import os
import rich.console
import statistics
from pathlib import Path
import Bio.Data.CodonTable

from Bio import SeqIO

from collections import OrderedDict, defaultdict

import taranis.utils


log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)


class AnalyzeSchema:
    def __init__(
        self,
        schema_allele: str,
        output: str,
        remove_subset: bool,
        remove_duplicated: bool,
        remove_no_cds: bool,
        genus: str,
        species: str,
        usegenus: str,
        prokka_cpus: int,
    ) -> "AnalyzeSchema":
        """AnalyzeSchema instance creation

        Args:
            self : AnalyzeSchema instance
            schema_allele (str): Folder path where schema files are located
            output (str): Out folder to save result
            remove_subset (bool): Remove subset sequences if True
            remove_duplicated (bool): Remove duplicated sequences if True
            remove_no_cds (bool): Removing non coding sequences if True
            genus (str): Genus name for Prokka schema genes annotation
            species (str): Species name for Prokka schema genes annotation
            usegenus (str): genus-specific BLAST databases for Prokka
            prokka_cpus (int): number of cpus used in prokka

        Returns:
            AnalyzeSchema: Instance of the created class
        """
        self.schema_allele = schema_allele
        self.allele_name = Path(self.schema_allele).stem
        self.output = output
        self.remove_subset = remove_subset
        self.remove_duplicated = remove_duplicated
        self.remove_no_cds = remove_no_cds
        self.genus = genus
        self.species = species
        self.usegenus = usegenus
        self.prokka_cpus = prokka_cpus

    def check_allele_quality(self, prokka_annotation: dict) -> OrderedDict:
        """Each allele in the locus file is analyzed its quality by checking
            if it can be converted to protein, has start/stop codon, has multiple
            stop codon, its a subsequence of the another allele, and if it is
            duplicated.
            A new schema file is generated, and if remove parameters are set, for
            the bad quality, they are removed in this new schema file.
            Dictionary with quality information for each allele is returned

        Args:
            self : AnalyzeSchema instance
            prokka_annotation (dict): Contains the annotation for each allele

        Returns:
            OrderedDict: Quality information for each allele
        """
        log.debug("Processing allele quality for %s", self.allele_name)
        a_quality = OrderedDict()
        allele_seq = {}
        bad_quality_record = []
        with open(self.schema_allele) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                try:
                    prokka_ann = prokka_annotation[record.id]["gene"]
                    product_annotation = prokka_annotation[record.id]["product"]
                except Exception:
                    prokka_ann = "Not found in prokka"
                    product_annotation = "Not found"
                a_quality[record.id] = {
                    "allele_name": self.allele_name,
                    "quality": "Good quality",
                    "reason": "-",
                    "direction": "forward",
                    "start_codon_alt": "standard",
                    "protein_seq": "",
                    "cds_coding": prokka_ann,
                    "product_annotation": product_annotation,                }
                allele_seq[record.id] = str(record.seq)
                a_quality[record.id]["length"] = str(len(str(record.seq)))
                a_quality[record.id]["dna_seq"] = str(record.seq)
                sequence_direction = taranis.utils.get_seq_direction(str(record.seq))

                if sequence_direction == "reverse":
                    record.seq = record.seq.reverse_complement()
                    a_quality[record.id]["direction"] = sequence_direction
                elif sequence_direction == "Error":
                    a_quality[record.id]["direction"] = "-"
                try:
                    a_quality[record.id]["protein_seq"] = str(
                        record.seq.translate(table=1, cds=True)
                    )

                except Bio.Data.CodonTable.TranslationError as e:
                    if "not a start codon" in str(e):
                        try:
                            # Check if sequence has an alternative start codon
                            # for protein coding
                            a_quality[record.id]["protein_seq"] = str(
                                record.seq.translate(table=2, cds=True)
                            )
                            a_quality[record.id]["start_codon_alt"] = "alternative"
                        except Bio.Data.CodonTable.TranslationError as e_2:
                            if "stop" in str(e_2):
                                a_quality[record.id]["reason"] = str(e_2).replace(
                                    "'", ""
                                )
                            else:
                                a_quality[record.id]["reason"] = str(e).replace("'", "")
                            a_quality[record.id]["quality"] = "Bad quality"
                    else:
                        a_quality[record.id]["quality"] = "Bad quality"
                        a_quality[record.id]["reason"] = str(e).replace("'", "")

                if (
                    self.remove_no_cds
                    and a_quality[record.id]["quality"] == "Bad quality"
                ):
                    bad_quality_record.append(record.id)
        log.debug("Checking bad quality of alleles for %s", self.allele_name)
        # check if there are duplicated alleles
        # get the unique sequences and compare the length with all sequences
        unique_seq = list(set(list(allele_seq.values())))
        value_to_keys = defaultdict(list)
        for rec_id, seq_value in allele_seq.items():
            value_to_keys[seq_value].append(rec_id)
            # Check if sequence is already duplicate
            if len(value_to_keys[seq_value]) > 1:
                a_quality[rec_id]["quality"] = "Bad quality"
                a_quality[rec_id]["reason"] = "Duplicate allele"
                if self.remove_duplicated:
                    bad_quality_record.append(rec_id)
            # check if sequence is a sub allele
            try:
                unique_seq.remove(seq_value)
            except ValueError:
                log.warning(
                    "Already  deleted same sequence as for record id  %s", record.id
                )
            if seq_value in unique_seq:
                a_quality[rec_id]["quality"] = "Bad quality"
                a_quality[rec_id]["reason"] = "Sub set allele"
                if self.remove_subset:
                    bad_quality_record.append(rec_id)

        new_schema_folder = os.path.join(self.output, "new_schema")
        _ = taranis.utils.create_new_folder(new_schema_folder)
        new_schema_file = os.path.join(new_schema_folder, self.allele_name + ".fasta")
        with open(self.schema_allele, "r") as _:
            with open(new_schema_file, "w") as fo:
                for record in SeqIO.parse(self.schema_allele, "fasta"):
                    if len(bad_quality_record) > 0:
                        if record.id not in bad_quality_record:
                            SeqIO.write(record, fo, "fasta")
                    else:
                        SeqIO.write(record, fo, "fasta")
        # update the schema allele with the new file
        self.schema_allele = new_schema_file

        return a_quality

    def fetch_statistics_from_alleles(self, a_quality: dict) -> dict:
        """By using the information for each allele in the input data create a
        dictionary with statistics data about length and quality

        Args:
            self: AnalyzeSchema instance
            a_quality (dict): Containing allele information

        Returns:
            dict: statistics information for all alleles
        """
        stderr.print("Processing quality statistics")
        record_data = {}
        bad_quality_reason = {}
        a_length = []
        bad_quality_counter = 0
        for record_id in a_quality.keys():
            record_data["allele_name"] = self.allele_name
            a_length.append(int(a_quality[record_id]["length"]))
            if a_quality[record_id]["quality"] == "Bad quality":
                bad_quality_counter += 1
                for reason in taranis.utils.POSIBLE_BAD_QUALITY:
                    if reason in a_quality[record_id]["reason"]:
                        bad_quality_reason[reason] = (
                            bad_quality_reason.get(reason, 0) + 1
                        )
        total_alleles = len(a_length)
        record_data["min_length"] = min(a_length)
        record_data["max_length"] = max(a_length)
        record_data["num_alleles"] = total_alleles
        record_data["mean_length"] = round(statistics.mean(a_length), 2)
        record_data["good_percent"] = round(
            100 * (total_alleles - bad_quality_counter) / total_alleles, 2
        )
        for item in taranis.utils.POSIBLE_BAD_QUALITY:
            record_data[item] = (
                bad_quality_reason[item] if item in bad_quality_reason else 0
            )

        return record_data

    def analyze_allele_in_schema(self) -> list[dict, dict]:
        """Analyze the alleles in the schema file by callig the function to
        annotate each of the alle and using this info to provide it for checking
        allele quality on check_allele_quality. With both info collection
        statistics are obtain, returning a list of 2 dict, one with the statistics
        information and the quality.

        Args:
            self : _description_

        Returns:
            list[dict, dict]: _description_
        """
        allele_data = {}
        log.info("Analizing allele %s", self.allele_name)
        # run annotations
        prokka_folder = os.path.join(self.output, "prokka", self.allele_name)
        anotation_files = taranis.utils.create_annotation_files(
            self.schema_allele, prokka_folder, self.allele_name, cpus=self.prokka_cpus
        )
        log.info("Fetching anotation information for %s", self.allele_name)
        prokka_annotation = taranis.utils.read_annotation_file(anotation_files + ".gff")

        # Perform quality
        a_quality = self.check_allele_quality(prokka_annotation)
        allele_data = self.fetch_statistics_from_alleles(a_quality)
        return [allele_data, a_quality]


def parallel_execution(
    schema_allele: str,
    output: str,
    remove_subset: bool,
    remove_duplicated: bool,
    remove_no_cds: bool,
    genus: str,
    species: str,
    usegenus: str,
    prokka_cpus: int = 3,
) -> list[dict, dict]:
    """_summary_

    Args:
        schema_allele (str): Folder path where schema files are located
        output (str): Out folder to save result
        remove_subset (bool): Remove subset sequences if True
        remove_duplicated (bool): Remove duplicated sequences if True
        remove_no_cds (bool): Removing non coding sequences if True
        genus (str): Genus name for Prokka schema genes annotation
        species (str): Species name for Prokka schema genes annotation
        usegenus (str): genus-specific BLAST databases for Prokka
        prokka_cpus (int): number of cpus used to execute prokka. Default 3

    Returns:
        list[dict, dict]:: _description_
    """
    schema_obj = AnalyzeSchema(
        schema_allele,
        output,
        remove_subset,
        remove_duplicated,
        remove_no_cds,
        genus,
        species,
        usegenus,
        prokka_cpus,
    )
    return schema_obj.analyze_allele_in_schema()


def collect_statistics(data, out_folder, output_allele_annot):
    def stats_graphics(stats_folder: str) -> None:
        """Create the statistics graphics. Pie graphic for allele quality,
            bar graphic for number of alleles, and box plot for allele variability

        Args:
            stats_folder (str): folder path to store graphic
        """
        stderr.print("Creating graphics")
        log.info("Creating graphics")
        allele_range = [0, 300, 600, 1000, 1500]
        graphic_folder = os.path.join(stats_folder, "graphics")
        _ = taranis.utils.create_new_folder(graphic_folder)

        # create graphic for alleles/number of genes
        group_alleles_df = stats_df.groupby(
            pd.cut(stats_df["num_alleles"], allele_range), observed=False
        ).count()
        _ = taranis.utils.create_graphic(
            graphic_folder,
            "num_genes_per_allele.png",
            "bar",
            allele_range[1:],
            group_alleles_df["num_alleles"].to_list(),
            ["Allele", "number of genes"],
            "Number of alleles per gene",
        )

        sum_all_alleles = stats_df["num_alleles"].sum()

        labels = taranis.utils.POSIBLE_BAD_QUALITY
        values = [stats_df[item].sum() for item in labels]

        labels.append("Good quality")
        values.append(sum_all_alleles - sum(values))
        _ = taranis.utils.create_graphic(
            graphic_folder,
            "quality_percent.png",
            "pie",
            values,
            "",
            labels,
            "Quality percent",
        )
        # create box plot for allele length variability
        _ = taranis.utils.create_graphic(
            graphic_folder,
            "allele_variability.png",
            "box",
            "",
            stats_df["mean_length"].to_list(),
            "",
            "Allele length variability",
        )

    summary_data = []
    a_quality = []
    for idx in range(len(data)):
        summary_data.append(data[idx][0])
        a_quality.append(data[idx][1])

    stats_df = pd.DataFrame(summary_data)
    stats_folder = os.path.join(out_folder, "statistics")
    _ = taranis.utils.create_new_folder(stats_folder)
    _ = taranis.utils.write_data_to_file(stats_folder, "statistics.csv", stats_df)
    stats_graphics(stats_folder)

    if output_allele_annot:
        # if parameter to save allele annotation then write to file
        ann_heading = [
            "gene",
            "allele",
            "allele direction",
            "nucleotide sequence",
            "protein sequence",
            "nucleotide sequence length",
            "star codon",
            "CDS coding",
            "product annotation",
            "allele quality",
            "bad quality reason",
        ]
        ann_fields = [
            "direction",
            "dna_seq",
            "protein_seq",
            "length",
            "start_codon_alt",
            "cds_coding",
            "product_annotation",
            "quality",
            "reason",
        ]

        ann_data = ",".join(ann_heading) + "\n"
        for gene in a_quality:
            for allele in gene.keys():
                data_field = [gene[allele][field] for field in ann_fields]
                ann_data += (
                    gene[allele]["allele_name"]
                    + ","
                    + allele
                    + ","
                    + ",".join(data_field)
                    + "\n"
                )
        _ = taranis.utils.write_data_to_compress_filed(
            out_folder, "allele_annotation.csv", ann_data
        )
    return
