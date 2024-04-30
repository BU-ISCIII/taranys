import io
import logging
import pandas as pd
import numpy as np
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
        self.allele_matrix = io.StringIO()
        self.allele_matrix.write("alleles\t" + "\t".join(allele_names) + "\n")
        self.allele_matrix.write("\n".join(out_data[1:]))
        self.allele_matrix.seek(0)
        matrix_pd = pd.read_csv(
            self.allele_matrix, sep="\t", index_col="alleles", engine="python"
        ).fillna(0)
        # Close object and discard memory buffer
        self.allele_matrix.close()
        log.debug(f"create distance for {allele_name}")
        return matrix_pd


class HammingDistance:
    def __init__(self, allele_matrix: pd.DataFrame) -> "HammingDistance":
        """HammingDistance instance creation

        Args:
            self.allele_matrix (pd.DataFrame): Distance matrix

        Returns:
            HammingDistance: created hamming distance
        """
        self.allele_matrix = allele_matrix

    def create_matrix(self, mask_values: list) -> pd.DataFrame:
        """Create hamming distance matrix

        Args:
            mask_values: list of values to mask p.e ["ASM", "LNF"]

        Returns:
            pd.DataFrame: Hamming distance matrix as panda DataFrame
        """
        # Mask unwanted values directly in the DataFrame
        regex_pattern = '|'.join([f".*{value}.*" for value in mask_values])
        self.allele_matrix.replace(regex_pattern, np.nan, regex=True, inplace=True)

        # Get unique values excluding NaN
        unique_values = pd.unique(
            self.allele_matrix.values.ravel("K")
        )
        unique_values = unique_values[~pd.isna(unique_values)]  # Exclude NaNs from unique values

        # Create binary matrix ('1' or '0' ) matching the input matrix vs the unique_values[0]
        # astype(int) is used to transform the boolean matrix into integer
        U = self.allele_matrix.eq(unique_values[0]).astype(int)
        # multiply the matrix with the transpose
        H = U.dot(U.T)

        # Repeat for each unique value
        for unique_val in range(1, len(unique_values)):
            U = self.allele_matrix.eq(unique_values[unique_val]).astype(int)
            # Add the value of the binary matrix with the previous stored values
            H = H.add(U.dot(U.T))

        # Convert to Boolean where True is not NaN (valid)
        valid_data = self.allele_matrix.notna()

        # Use broadcasting to find pairwise non-NaN entries
        # valid_data[:, None] adds a new axis, making it a 3D array where each 2D slice is one sample's valid data
        # We then logical AND across all pairs of samples
        pairwise_valid = valid_data.values[:, None] & valid_data.values

        # Sum along the third dimension to get pairwise counts of non-NaN positions
        pairwise_valid_counts = pairwise_valid.sum(axis=2)
        distance_matrix = pairwise_valid_counts - H

        return pd.DataFrame(distance_matrix, index=self.allele_matrix.index, columns=self.allele_matrix.index)
