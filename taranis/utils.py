#!/usr/bin/env python
import Bio.Data.CodonTable
import glob
import gzip
import io
import logging
import multiprocessing
import os
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import questionary
import shutil

import re
import rich.console

import subprocess
import tarfile

import sys

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict

import warnings
from Bio import BiopythonWarning

# import pdb

log = logging.getLogger(__name__)


def rich_force_colors():
    """
    Check if any environment variables are set to force Rich to use coloured output
    """
    if (
        os.getenv("GITHUB_ACTIONS")
        or os.getenv("FORCE_COLOR")
        or os.getenv("PY_COLORS")
    ):
        return True
    return None


stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=rich_force_colors(),
)

START_CODON_FORWARD = ["ATG", "ATA", "ATT", "GTG", "TTG"]
START_CODON_REVERSE = ["CAT", "TAT", "AAT", "CAC", "CAA"]

STOP_CODON_FORWARD = ["TAA", "TAG", "TGA"]
STOP_CODON_REVERSE = ["TTA", "CTA", "TCA"]

POSIBLE_BAD_QUALITY = [
    "not a start codon",
    "not a stop codon",
    "Extra in frame stop codon",
    "is not a multiple of three",
    "Duplicate allele",
    "Sub set allele",
]


def has_start_codon(seq):
    return seq[:3] in START_CODON_FORWARD or seq[-3:] in START_CODON_REVERSE


def has_stop_codon(seq):
    return seq[:3] in STOP_CODON_FORWARD or seq[-3:] in STOP_CODON_REVERSE


def cpus_available() -> int:
    """Get the number of cpus available in the system

    Returns:
        int: number of cpus
    """
    return multiprocessing.cpu_count()


def get_seq_direction(allele_sequence):
    # check direction
    if (
        allele_sequence[0:3] in START_CODON_FORWARD
        or allele_sequence[-3:] in STOP_CODON_FORWARD
    ):
        return "forward"
    if (
        allele_sequence[-3:] in START_CODON_REVERSE
        or allele_sequence[0:3] in STOP_CODON_REVERSE
    ):
        return "reverse"
    return "both"


def check_additional_programs_installed(software_list: list) -> None:
    """Check if the input list of programs are installed in the system

    Args:
        software_list (list): list of programs to be checked
    """
    for program, command in software_list:
        try:
            _ = subprocess.run(
                [program, command],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
        except Exception as e:
            log.error(
                "Program %s is not installed in the system. Error message: %s ",
                program,
                e,
            )
            stderr.print("[red] Program " + program + " is not installed in the system")
            sys.exit(1)
    return


def convert_to_protein(sequence: str, force_coding: bool = False) -> dict:
    """Check if the input sequence is a coding protein.

    Args:
        sequence (str): sequence to be checked
        force_coding (bool, optional): force to check if sequence is coding.
            Defaults to False.

    Returns:
        direction(str): reverse or forward
        protein (str): protein sequence
        error (bool): True/False
        error_detail (str): translate method error
    """
    protein = "-"
    error = False
    error_detail = "-"

    direction = get_seq_direction(sequence)

    seq = Seq(sequence)

    if direction == "reverse":
        seq = seq.reverse_complement()
    try:
        # Table 11 is for bacteria, archaea and chloroplast
        protein = seq.translate(table=11, to_stop=False, cds=force_coding)
    except Bio.Data.CodonTable.TranslationError as e:
        error = True
        error_detail = str(e)
        log.debug(f"Error when translating protein {error_detail}")

    return direction, str(protein), error, error_detail


def create_annotation_files(
    fasta_file: str,
    annotation_dir: str,
    prefix: str,
    genus: str = "Genus",
    species: str = "species",
    usegenus: str = False,
    cpus: int = 3,
) -> str:
    """prokka command is executed to generate the annotation files.
    Return the folder path where prokka store these files

    Args:
        fasta_file (str): fasta file used for annotation
        annotation_dir (str): folder where annotation files are saved
        prefix (str): string used for naming annotation files
        genus (str, optional): parameter used in proka. Defaults to "Genus".
        species (str, optional): parameter used in proka. Defaults to "species".
        usegenus (str, optional): _description_. Defaults to False.
        cpus (int, optional): number of cpus used to run prokka. Defaults to 3.

    Returns:
        str: folder path where generated files from prokka are stored
    """
    try:
        _ = subprocess.run(
            [
                "prokka",
                fasta_file,
                "--force",
                "--outdir",
                annotation_dir,
                "--genus",
                genus,
                "--species",
                species,
                "--usegenus",
                str(usegenus),
                "--gcode",
                "11",
                "--prefix",
                prefix,
                "--cpus",
                str(cpus),
                "--quiet",
            ]
        )
    except Exception as e:
        log.error("Unable to run prokka. Error message: %s ", e)
        stderr.print("[red] Unable to run prokka. Given error; " + e)
        sys.exit(1)
    return os.path.join(annotation_dir, prefix)


def create_new_folder(folder_name: str) -> None:
    """Create directory defined in input data. No error occurs if folder exists

    Args:
        folder_name (str): folder path to be created
    """
    try:
        os.makedirs(folder_name, exist_ok=True)
    except Exception as e:
        log.error("Folder %s can not be created %s", folder_name, e)
        stderr.print("[red] Folder does not have any file which match your request")
        sys.exit(1)
    return


def create_graphic(
    out_folder: str,
    f_name: str,
    mode: str,
    x_data: list,
    y_data: list,
    labels: list,
    title: str,
) -> None:
    """Create the graphic and save it to file

    Args:
        out_folder (str): folder path to save the graphic
        f_name (str): file name including extension
        mode (str): type of graphic
        x_data (list): data for x axis
        y_data (list): data for y axis
        labels (list): labels to be included
        title (str): title of the figure
    """
    fig = go.Figure()
    layout_update = {}
    plot_options = {
        "lines": (go.Scatter, {"mode": mode}),
        "pie": (go.Pie, {"labels": labels, "values": x_data}),
        "bar": (go.Bar, {"x": x_data, "y": y_data}),
        "box": (go.Box, {"y": y_data}),
    }

    if mode in plot_options:
        trace_class, trace_kwargs = plot_options[mode]
        fig.add_trace(trace_class(**trace_kwargs))
        if mode == "bar":
            layout_update = {
                "xaxis_title": labels[0],
                "yaxis_title": labels[1],
                "xaxis_tickangle": 45,
            }
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    layout_update["title_text"] = title
    fig.update_layout(**layout_update)

    pio.write_image(fig, os.path.join(out_folder, f_name))
    return


def delete_folder(folder_to_delete: str) -> None:
    """Delete the input folder

    Args:
        folder_to_delete (str): folder path to be deleted
    """
    try:
        shutil.rmtree(folder_to_delete)
    except Exception as e:
        log.error("Folder %s can not be deleted %s", folder_to_delete, e)
        stderr.print("[red] Folder does not have any file which match your request")
        sys.exit(1)
    return


def file_exists(file_to_check):
    """Checks if input file exists

    Args:
        file_to_check (string): file name  including path of the file

    Returns:
        boolean: True if exists
    """
    if os.path.isfile(file_to_check):
        return True
    return False


"""
def find_nearest_numpy_value(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

 """


def filter_data_frame_by_parameters(
    data_frame: pd.DataFrame,
    column_thr: int,
    row_thr: int,
    filter_str: list[str],
    replaced_by_zero: bool,
) -> pd.DataFrame:
    # get the number of columns and rows
    num_rows, num_columns = data_frame.shape
    # remove the columns which the filter strings are higher than the threshold
    column_threshold = column_thr * num_rows / 100
    # Condition: Check if any string in the filter list is present in each cell of the DataFrame
    f_condition = data_frame.apply(
        lambda column: column.astype(str).str.contains("|".join(filter_str), na=False)
    )
    if replaced_by_zero:
        new_data_frame = data_frame.mask(f_condition, 0)
    else:
        # Count the number of hits per column
        hits_per_column = f_condition.sum()
        # pdb.set_trace()
        # Filter for removing columns where the count of hits is higher than the threshold
        to_be_removed_columns = hits_per_column[
            hits_per_column > column_threshold
        ].index
        new_data_frame = data_frame.drop(columns=to_be_removed_columns)
    # pdb.set_trace()
    # remove the rows which the filter strings are higher than the threshold
    # row_threshold = row_thr * num_columns
    return new_data_frame


def folder_exists(folder_to_check):
    """Checks if input folder exists

    Args:
        folder_to_check (string): folder name  including path

    Returns:
        boolean: True if exists
    """
    if os.path.isdir(folder_to_check):
        return True
    return False


def get_alignment_data(ref_sequence: str, allele_sequence: str, ref_allele) -> dict:
    """Get the alignment data between the reference allele and the match allele
        sequence. It returns 3 lines, the reference allele, the alignment character
        and the match allele sequence

    Args:
        allele_sequence (str): sequence to be compared
        ref_sequences (dict): sequences of reference alleles

    Returns:
        dict: key: ref_sequence, value: alignment data
    """
    alignment_data = {}
    alignment = ""
    for _, (ref, alt) in enumerate(zip(ref_sequence, allele_sequence)):
        if ref == alt:
            alignment += "|"
        else:
            alignment += " "
    alignment_data[ref_allele] = [ref_sequence, alignment, allele_sequence]
    return alignment_data


def get_files_in_folder(folder: str, extension: str = None) -> list[str]:
    """get the list of files, filtered by extension in the input folder. If
    extension is not set, then all files in folder are returned

    Args:
        folder (str): Folder path
        extension (str, optional): Extension for filtering. Defaults to None.

    Returns:
        list[str]: list of files which match the condition
    """
    if not folder_exists(folder):
        log.error("Folder %s does not exists", folder)
        stderr.print("[red] Schema folder does not exist. " + folder + "!")
        sys.exit(1)
    if extension is None:
        extension = "*"
    folder_files = os.path.join(folder, "*." + extension)
    return glob.glob(folder_files)


def get_multiple_alignment(input_buffer: io.StringIO, mafft_cpus: int) -> list[str]:
    """Run MAFFT with input from the string buffer and capture output to another string buffer

    Args:
        input_buffer (io.StringIO): fasta sequences to be aligned
        mafft_cpus (int): number of cpus to be used in mafft
    Returns:
        list[str]: list of aligned sequences
    """
    output_buffer = io.StringIO()
    # Run MAFFT
    mafft_command = (
        "mafft --auto --quiet --thread " + str(mafft_cpus) + " -"
    )  # "-" tells MAFFT to read from stdin
    process = subprocess.Popen(
        mafft_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE
    )
    stdout, _ = process.communicate(input_buffer.getvalue().encode())

    # Convert the stdout bytes to a string buffer
    output_buffer = io.StringIO(stdout.decode())
    output_buffer.seek(0)
    # convert the string buffer to a list of lines
    multi_result = []
    for line in output_buffer:
        multi_result.append(line)

    # Close the file objects and process
    output_buffer.close()
    process.stdout.close()

    return multi_result


def get_snp_information(
    ref_sequence: str, alt_sequence: str, ref_allele_name
) -> dict[list[str]]:
    """Get the snp information between the reference allele sequence and the
        allele sequence in sample.
        It collects; position of snp, nucleotide changed reference/alternative,
        triplet code (belongs the change), amino acid change and category of
        amino acid

    Args:
        ref_sequences (dict): sequences of reference alleles
        allele_sequence (str): sequence to be compared

    Returns:
        dict: key: ref_sequence, value: list of snp information
    """
    # Supress warning that len of alt sequence  not a multiple of three
    warnings.simplefilter("ignore", BiopythonWarning)
    snp_info = {}
    ref_protein = str(Seq(ref_sequence).translate())
    if len(alt_sequence) % 3 != 0:
        log.debug(
            "Sequence %s is not a multiple of three. Removing last nucleotides",
            ref_allele_name,
        )
        # remove the last nucleotides to be able to translate to protein
        alt_sequence = alt_sequence[: len(alt_sequence) // 3 * 3]

    alt_protein = str(Seq(alt_sequence).translate())
    snp_line = []
    # get the shortest sequence for the loop
    length_for_snp = min(len(ref_sequence), len(alt_sequence))
    for idx in range(length_for_snp):
        if ref_sequence[idx] != alt_sequence[idx]:
            # calculate the triplet index
            triplet_idx = idx // 3
            # get triplet code
            ref_triplet = ref_sequence[triplet_idx * 3 : triplet_idx * 3 + 3]
            alt_triplet = alt_sequence[triplet_idx * 3 : triplet_idx * 3 + 3]
            # get amino acid change
            ref_aa = ref_protein[triplet_idx]
            try:
                alt_aa = alt_protein[triplet_idx]
            except IndexError as e:
                log.debug(
                    "Unable to get amino acid for %s and %s with error %s",
                    ref_allele_name,
                    alt_sequence,
                    e,
                )
                alt_aa = "-"
            # get amino acid category
            ref_category = map_amino_acid_to_annotation(ref_sequence[triplet_idx])
            alt_category = map_amino_acid_to_annotation(alt_sequence[triplet_idx])
            snp_line.append(
                [
                    str(idx),
                    ref_sequence[idx],
                    alt_sequence[idx],
                    ref_triplet,
                    alt_triplet,
                    ref_aa,
                    alt_aa,
                    ref_category,
                    alt_category,
                ]
            )
    snp_info[ref_allele_name] = snp_line
    return snp_info


def map_amino_acid_to_annotation(amino_acid):
    # Dictionary mapping amino acids to their categories
    amino_acid_categories = {
        "A": "Nonpolar",
        "C": "Polar",
        "D": "Acidic",
        "E": "Acidic",
        "F": "Nonpolar",
        "G": "Nonpolar",
        "H": "Basic",
        "I": "Nonpolar",
        "K": "Basic",
        "L": "Nonpolar",
        "M": "Nonpolar",
        "N": "Polar",
        "P": "Nonpolar",
        "Q": "Polar",
        "R": "Basic",
        "S": "Polar",
        "T": "Polar",
        "V": "Nonpolar",
        "W": "Nonpolar",
        "Y": "Polar",
    }

    # Return the category of the given amino acid
    return amino_acid_categories.get(amino_acid, "Unknown")


def prompt_text(msg):
    source = questionary.text(msg).unsafe_ask()
    return source


def prompt_user_if_folder_exists(folder: str) -> bool:
    """Prompt the user to continue if the folder exists

    Args:
        folder (str): folder path

    Returns:
        bool: True if user wants to continue
    """
    if folder_exists(folder):
        q_question = (
            "Folder "
            + folder
            + " already exists. Files will be overwritten. Do you want to continue?"
        )
        if "no" in query_user_yes_no(q_question, "no"):
            log.info("Aborting code by user request")
            stderr.print("[red] Exiting code. ")
            sys.exit(1)
    else:
        try:
            os.makedirs(folder)
        except OSError as e:
            log.info("Unable to create folder at %s with error %s", folder, e)
            stderr.print("[red] ERROR. Unable to create folder  " + folder)
            sys.exit(1)

    return True


def query_user_yes_no(question, default):
    """Query the user to choose yes or no for the query question

    Args:
        question (string): Text message
        default (string): default option to be used: yes or no

    Returns:
        user select: True continue with code
    """
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)
    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == "":
            return valid[default]
        elif choice in valid:
            if "y" in choice:
                return "yes"
            else:
                return "no"
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")


def read_annotation_file(ann_file: str) -> dict:
    """Read the annotation file and return a dictionary where key is the allele
    name and the value is the annotation data that prokka was defined for the
    allele

    Args:
        ann_file (str): annotation file path (gff)

    Returns:
        dict: contains the allele name and the predction
    """
    """example of annotation file

    lmo0002_782	Prodigal:002006	CDS	1	1146	.	+	0	ID=OJGEGONH_00782;Name=dnaN_782;db_xref=COG:COG0592;gene=dnaN_782;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P05649;locus_tag=OJGEGONH_00782;product=Beta sliding clamp
    lmo0002_783	Prodigal:002006	CDS	1	1146	.	+	0	ID=OJGEGONH_00783;Name=dnaN_783;db_xref=COG:COG0592;gene=dnaN_783;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P05649;locus_tag=OJGEGONH_00783;product=Beta sliding clamp
    lmo0049_3	Prodigal:002006	CDS	1	162	.	+	0	ID=CODOCEEL_00001;inference=ab initio prediction:Prodigal:002006;locus_tag=CODOCEEL_00001;product=hypothetical protein
    lmo0049_6	Prodigal:002006	CDS	1	162	.	+	0	ID=CODOCEEL_00002;inference=ab initio prediction:Prodigal:002006;locus_tag=CODOCEEL_00002;product=hypothetical protein

    """
    ann_data = {}
    with open(ann_file, "r") as fh:
        lines = fh.readlines()

    for line in lines:
        if "Prodigal" in line:
            gene_match = re.search(r"(.*)[\t]Prodigal.*gene=(\w+)_.*product=(.*)", line)
            if gene_match:
                ann_data[gene_match.group(1)] = {
                    "gene": gene_match.group(2),
                    "product": gene_match.group(3).strip(),
                }
            else:
                pred_match = re.search(r"(.*)[\t]Prodigal.*product=(\w+)_.*", line)
                if pred_match:
                    ann_data[pred_match.group(1)] = pred_match.group(2).strip()
        if "fasta" in line:
            break
    return ann_data


def read_compressed_file(
    file_name: str, separator: str = ",", index_key: int = None, mapping: list = []
) -> dict | str:
    """Read the compressed file and return a dictionary using as key value
    the mapping data if the index_key is an integer, else return the uncompressed
    file

    Args:
        file_name (str): file to be uncompressed
        separator (str, optional): split line according separator. Defaults to ",".
        index_key (int, optional): index value . Defaults to None.
        mapping (list, optional): defined the key value for dictionary. Defaults to [].

    Returns:
        dict|str: uncompresed information file
    """
    out_data = {}
    with gzip.open(file_name, "rb") as fh:
        lines = fh.readlines()
    if not index_key:
        return lines[:-2]
    for line in lines[1:]:
        line = line.decode("utf-8")
        s_line = line.split(separator)
        # ignore empty lines
        if len(s_line) == 1:
            continue
        key_data = s_line[index_key]
        out_data[key_data] = {}
        for item in mapping:
            out_data[key_data][item[0]] = s_line[item[1]]
    return out_data


def read_fasta_file(fasta_file: str, convert_to_dict=False) -> dict | str:
    """Read the fasta file and return the data as a dictionary if convert_to_dict

    Args:
        fasta_file (str): _description_
        convert_to_dict (bool, optional): _description_. Defaults to False.

    Returns:
        dict: fasta id as key and sequence as value in str format
    """
    conv_fasta = OrderedDict()
    if convert_to_dict:
        with open(fasta_file, "r") as fh:
            for record in SeqIO.parse(fh, "fasta"):
                conv_fasta[record.id] = str(record.seq)
        return conv_fasta
    return SeqIO.parse(fasta_file, "fasta")


def write_fasta_file(
    out_folder: str, f_name: str, allele_name: str, seq_data: str
) -> str:
    """_summary_

    Args:
        out_folder (str): _description_
        seq_data (str): _description_
        allele_name (str, optional): _description_. Defaults to None.
        f_name (str, optional): _description_. Defaults to None.

    Returns:
        str: _description_
    """
    try:
        os.makedirs(out_folder, exist_ok=True)
    except OSError as e:
        print(e)
        sys.exit(1)

    f_path_name = os.path.join(out_folder, f_name + ".fasta")
    with open(f_path_name, "w") as fo:
        fo.write("> " + allele_name + "\n")
        fo.write(seq_data)
    return f_path_name


def write_data_to_compress_filed(out_folder, f_name, dump_data):
    with io.BytesIO() as buffer:
        with tarfile.open(fileobj=buffer, mode="w:gz") as tar:
            # Add data to the tar archive
            tarinfo = tarfile.TarInfo(f_name)
            # Example: Write a string to the tar.gz file (replace this with your data)
            data_bytes = dump_data.encode("utf-8")
            tarinfo.size = len(data_bytes)
            tar.addfile(tarinfo, io.BytesIO(data_bytes))

        # Get the content of the in-memory tar.gz file
        buffer.seek(0)
        tar_data = buffer.read()
    file_path_name = os.path.join(out_folder, Path(f_name).stem + ".tar.gz")
    with open(file_path_name, "wb") as fo:
        fo.write(tar_data)


def write_data_to_file(
    out_folder: str,
    f_name: str,
    data: pd.DataFrame | list,
    include_header: bool = True,
    data_type: str = "pandas",
    extension: str = "csv",
) -> None:
    """write data in the input parameter to disk

    Args:
        out_folder (str): Folder path to store file
        f_name (str): file name without extension
        data (pd.DataFrame | list): data to write. Can be dataframe or list
        include_header (bool, optional): for pandas input check if header has to
            be included in file. Defaults to True.
        data_type (str, optional): type of data pandas or list. Defaults to "pandas".
        extension (str, optional): extension of file. Defaults to "csv".
    """
    f_path_name = os.path.join(out_folder, f_name)
    if data_type == "pandas":
        data.to_csv(f_path_name, sep=",", header=include_header)
        return


"""
def find_multiple_stop_codons(seq) :
    stop_codons = ['TAA', 'TAG','TGA']
    c_index = []
    for idx in range (0, len(seq) -2, 3) :
        c_seq = seq[idx:idx + 3]
        if c_seq in stop_codons :
            c_index.append(idx)
    if len(c_index) == 1:
        return False
    return True
"""
