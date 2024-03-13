#!/usr/bin/env python
import glob
import gzip
import io
import logging
import multiprocessing
import numpy as np
import os
import pandas as pd
import plotly.graph_objects as go
import questionary
import shutil

import re
import rich.console

import subprocess
import tarfile

import sys

from pathlib import Path
from Bio import SeqIO

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


def cpus_available() -> int:
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
    return "Error"


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
    if mode == "lines":
        fig.add_trace(go.Scatter(x=x_data, y=y_data, mode=mode, name=title))
        fig.update_layout(xaxis_title=labels[0], yaxis_title=labels[1])
    elif mode == "pie":
        fig.add_trace(go.Pie(labels=labels, values=x_data))
    elif mode == "bar":
        fig.add_trace(go.Bar(x=x_data, y=y_data))
        fig.update_layout(xaxis_title=labels[0], yaxis_title=labels[1])
    elif mode == "box":
        fig.add_trace(go.Box(y=y_data))

    fig.update_layout(title_text=title)
    fig.write_image(os.path.join(out_folder, f_name))
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


def find_nearest_numpy_value(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


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


def grep_execution(input_file: str, pattern: str, parameters: str) -> list:
    """_summary_

    Args:
        input_file (str): _description_
        pattern (str): _description_
        parmeters (str): _description_

    Returns:
        list: _description_
    """
    try:
        result = subprocess.run(
            ["grep", parameters, pattern, input_file],
            capture_output=True,
            check=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        log.error("Unable to run grep. Error message: %s ", e)
        return []
    return result.stdout.split("\n")


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
                ann_data[gene_match.group(1)] = {"gene": gene_match.group(2) , "product": gene_match.group(3).strip()}
            else:
                pred_match = re.search(r"(.*)[\t]Prodigal.*product=(\w+)_.*", line)
                if pred_match:
                    ann_data[pred_match.group(1)] = pred_match.group(2).strip()
        if "fasta" in line:
            break
    return ann_data

def read_compressed_file(file_name: str, separator: str = ",", index_key: int=None, mapping: list=[]) -> dict|str:
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
        key_data =s_line[index_key]
        out_data[key_data] = {}
        for item in mapping:
            out_data[key_data][item[0]] = s_line[item[1]]
    return out_data

def read_fasta_file(fasta_file):
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
