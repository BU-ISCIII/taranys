import io
import logging
import pandas as pd
import subprocess
import rich
import sys
from pathlib import Path
import taranis.utils

log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)


class DistanceMatrix:
    def __init__(
        self, file_path: str, k_mer_value: str = "17", sketch_size: str = "2000"
    ) -> "DistanceMatrix":
        """DistanceMatrix instance creation

        Args:
            file_path (str): Locus file path
            k_mer_value (str, optional): Hashes will be based on strings of this many nucleotides. Defaults to "21".
            sketch_size (str, optional): Each sketch will have at most this many non-redundant min-hashes. Defaults to "2000".

        Returns:
            DistanceMatrix: created distance
        """
        self.file_path = file_path
        self.k_mer_value = k_mer_value
        self.sketch_size = sketch_size

    def create_matrix(self) -> pd.DataFrame:
        """Create distance matrix using external program called mash

        Returns:
            pd.DataFrame: Triangular distance matrix as panda DataFrame
        """
        allele_name = Path(self.file_path).stem
        mash_distance_command = [
            "mash",
            "triangle",
            "-i",
            self.file_path,
            "-k",
            str(self.k_mer_value),
            "-s",
            str(self.sketch_size),
        ]
        try:
            mash_distance_result = subprocess.Popen(
                mash_distance_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            out, _ = mash_distance_result.communicate()
            log.debug(f"calculate mash distance for {allele_name}")
        except Exception as e:
            log.error(f"Unable to create distance matrix for {self.file_path}. {e}")
            stderr.print(
                f"[red] Error when creating distance matrix for {self.file_path}"
            )
            stderr.print(f"{e}")
            sys.exit(1)
        finally:
            # Close the file handles
            mash_distance_result.stdout.close()
            mash_distance_result.stderr.close()

        out_data = out.decode("UTF-8").split("\n")
        allele_names = [item.split("\t")[0] for item in out_data[1:-1]]
        # create file in memory to increase speed
        dist_matrix = io.StringIO()
        dist_matrix.write("alleles\t" + "\t".join(allele_names) + "\n")
        dist_matrix.write("\n".join(out_data[1:]))
        dist_matrix.seek(0)
        matrix_pd = pd.read_csv(
            dist_matrix, sep="\t", index_col="alleles", engine="python"
        ).fillna(0)
        # Close object and discard memory buffer
        dist_matrix.close()
        log.debug(f"create distance for {allele_name}")
        return matrix_pd


class HammingDistance:
    def __init__(self, dist_matrix: pd.DataFrame) -> "HammingDistance":
        """HammingDistance instance creation

        Args:
            dist_matrix (pd.DataFrame): Distance matrix

        Returns:
            HammingDistance: created hamming distance
        """
        self.dist_matrix = dist_matrix

    def create_matrix(self) -> pd.DataFrame:
        """Create hamming distance matrix using external program called mash

        Returns:
            pd.DataFrame: Hamming distance matrix as panda DataFrame
        """

        unique_values = pd.unique(
            self.dist_matrix[list(self.dist_matrix.keys())].values.ravel("K")
        )
        # Create binary matrix ('1' or '0' ) matching the input matrix vs the unique_values[0]
        # astype(int) is used to transform the boolean matrix into integer
        U = self.dist_matrix.eq(unique_values[0]).astype(int)
        # multiply the matrix with the transpose
        H = U.dot(U.T)

        # Repeat for each unique value
        for unique_val in range(1, len(unique_values)):
            U = self.dist_matrix.eq(unique_values[unique_val]).astype(int)
            # Add the value of the binary matrix with the previous stored values
            H = H.add(U.dot(U.T))

        return len(self.dist_matrix.columns) - H

        """
         dist_matrix = self.dist_matrix
        allele_names = dist_matrix.index
        hamming_matrix = pd.DataFrame(index=allele_names, columns=allele_names)
        for i in allele_names:
            for j in allele_names:
                hamming_matrix.at[i, j] = sum(
                    dist_matrix.loc[i] != dist_matrix.loc[j]
                )
        return hamming_matrix
        """
