import logging
import os
import rich
import subprocess
import taranys.utils

from pathlib import Path
from Bio.Blast.Applications import NcbiblastnCommandline

log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranys.utils.rich_force_colors(),
)


class Blast:
    def __init__(self, db_type: str):
        """Blast instance creation

        Args:
            db_type (str): type of blast database (nucleotide or protein)
        """
        self.db_type = db_type

    def create_blastdb(self, file_name: str, blast_dir: str) -> None:
        """Create blast database and store it at blast dir

        Args:
            file_name (str): Fasta file from generate the database
            blast_dir (str): directory to store blast database files
        """
        self.f_name = Path(file_name).stem
        db_dir = os.path.join(blast_dir, self.f_name)
        self.out_blast_dir = os.path.join(db_dir, self.f_name)

        blast_command = [
            "makeblastdb",
            "-in",
            file_name,
            "-parse_seqids",
            "-dbtype",
            self.db_type,
            "-out",
            self.out_blast_dir,
        ]
        try:
            _ = subprocess.run(
                blast_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
        except Exception as e:
            log.error("Unable to create blast db for %s ", self.f_name)
            log.error(e)
            stderr.print(
                f"[red] Unable to create blast database for sample {self.f_name}"
            )
            exit(1)
        return

    def run_blast(
        self,
        query: str,
        evalue: float = 10,
        perc_identity: int = 90,
        reward: int = 1,
        penalty: int = -2,
        gapopen: int = 1,
        gapextend: int = 1,
        max_target_seqs: int = 10000,
        max_hsps: int = 1,
        num_threads: int = 1,
        query_type: str = "file",
    ) -> list:
        """blast command is executed, returning a list of each match found

        Args:
            query (str): file path which contains the fasta sequence to query
            evalue (float, optional): filtering results on e-value. Defaults to 0.001.
            perc_identity (int, optional): percentage of identity. Defaults to 90.
            reward (int, optional): value for rewardin. Defaults to 1.
            penalty (int, optional): penalty value. Defaults to -2.
            gapopen (int, optional): value for gap open. Defaults to 1.
            gapextend (int, optional): value for gap extended. Defaults to 1.
            max_target_seqs (int, optional): max target to output. Defaults to 1000.
            max_hsps (int, optional): max hsps. Defaults to 10.
            num_threads (int, optional): number of threads. Defaults to 1.
            query_type (str, optional): format of query (either file or string)
        Returns:
            list: list of strings containing blast results
        """
        if query_type == "stdin":
            stdin_query = query
            query = "-"
        blast_parameters = '"6 , qseqid , sseqid , pident ,  qlen , length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend , sseq , qseq , slen"'
        cline = NcbiblastnCommandline(
            task="blastn",
            db=self.out_blast_dir,
            evalue=evalue,
            perc_identity=perc_identity,
            reward=reward,
            penalty=penalty,
            gapopen=gapopen,
            gapextend=gapextend,
            outfmt=blast_parameters,
            max_target_seqs=max_target_seqs,
            max_hsps=max_hsps,
            num_threads=num_threads,
            query=query,
        )
        try:
            if query_type == "stdin":
                out, _ = cline(stdin=stdin_query)
            else:
                out, _ = cline()
        except Exception as e:
            log.error("Unable to run blast for %s ", self.out_blast_dir)
            log.error(e)
            stderr.print(f"[red] Unable to run blast {self.out_blast_dir}")
            exit(1)
        return out.splitlines()
