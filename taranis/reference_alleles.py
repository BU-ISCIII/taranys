import logging
import numpy as np
import rich.console
from pathlib import Path
import os

import taranis.utils
import taranis.distance
import taranis.clustering
import taranis.eval_cluster
from Bio import SeqIO

log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)


class ReferenceAlleles:
    def __init__(
        self,
        fasta_file: str,
        output: str,
        eval_cluster: bool,
        kmer_size: int,
        sketch_size: int,
        cluster_resolution: float = 0.75,
        seed: int = None,
    ):
        """ReferenceAlleles instance creation

        Args:
            fasta_file (str): file name included path for locus
            output (str): output folder
            eval_cluster (bool): True if cluster evaluation must be done
            kmer_size (int): kmer size for mash distance
            sketch_size (int): sketch size for mash distance
            cluster_resolution (float): resolution for clustering
            seed (int): seed for random number generator
        """
        self.fasta_file = fasta_file
        self.locus_name = Path(fasta_file).stem
        self.output = output
        self.eval_cluster = eval_cluster
        self.kmer_size = kmer_size
        self.sketch_size = sketch_size
        self.cluster_resolution = cluster_resolution
        self.seed = seed
        self.selected_locus = {}
        self.cluster_obj = None

    def create_distance_matrix(self) -> list:
        """Create the distance matrix for the alleles in the fasta file

        Returns:
            np.array: distance matrix
            dict: position to allele name
        """
        log.debug("Processing distance matrix for $s", self.fasta_file)
        distance_obj = taranis.distance.DistanceMatrix(
            self.fasta_file, self.kmer_size, self.sketch_size
        )
        mash_distance_df = distance_obj.create_matrix()
        log.debug(f"Created distance matrix for {self.fasta_file}")
        position_to_allele = {
            x: mash_distance_df.columns[x] for x in range(len(mash_distance_df.columns))
        }
        # convert the  triangle matrix into full data matrix
        matrix_np = mash_distance_df.to_numpy()
        t_matrix_np = matrix_np.transpose()
        matrix_np = t_matrix_np + matrix_np
        # At this point minimal distance is 0. For clustering requires to be 1
        # the oposite.
        dist_matrix_np = (matrix_np - 1) * -1
        return dist_matrix_np, position_to_allele

    def processing_cluster_data(
        self, cluster_data: np.array, cluster_ptrs: np.array, position_to_allele: dict
    ) -> dict:
        reference_alleles = []
        for cluster_id, values in cluster_data.items():
            center_allele = position_to_allele[values["center_id"]]
            values["center_id"] = center_allele
            reference_alleles.append(center_allele)
        alleles_in_cluster = self.cluster_obj.convert_to_seq_clusters(
            cluster_ptrs, position_to_allele
        )
        cluster_folder = os.path.join(self.output, "Clusters")
        _ = taranis.utils.create_new_folder(cluster_folder)
        cluster_file = os.path.join(
            cluster_folder, "cluster_alleles_" + self.locus_name + ".txt"
        )
        with open(cluster_file, "w") as fo:
            for cluster_id, alleles in alleles_in_cluster.items():
                fo.write("Cluster number" + str(cluster_id + 1) + "\n")
                fo.write("\n".join(alleles) + "\n")

        return {
            "cluster_data": cluster_data,
            "reference_alleles": reference_alleles,
            "alleles_in_cluster": alleles_in_cluster,
        }

    def save_reference_alleles(self, reference_alleles: list) -> str:
        """From the input list it fetch the allele squence and save it as fasta

        Args:
            reference_alleles (list): list having the allele ids

        Returns:
            str: file path of the reference alleles
        """
        record_seq = {}
        with open(self.fasta_file) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                if record.id in reference_alleles:
                    record_seq[record.id] = str(record.seq)
        ref_allele_file = os.path.join(self.output, self.locus_name + ".fasta")
        with open(ref_allele_file, "w") as fo:
            for ref_allele in reference_alleles:
                fo.write(">" + ref_allele + "\n")
                fo.write(record_seq[ref_allele] + "\n")
        return ref_allele_file

    def create_ref_alleles(self) -> dict:
        """Alleles in fasta file are clustering by using two additional classes:
            DistanceMatrix which creates a matrix of distance using the allele
            sequences, and ClusterDistance which get the matrix and group the
            alleles in clusters.

        Returns:
            dict: containg statistics information for each cluster, and
            optionally a list of evaluation cluster results
        """
        self.records = taranis.utils.read_fasta_file(self.fasta_file)
        dist_matrix_np, position_to_allele = self.create_distance_matrix()
        self.cluster_obj = taranis.clustering.ClusterDistance(
            dist_matrix_np,
            self.locus_name,
        )

        for resolution in np.arange(self.cluster_resolution, 1, 0.025):
            cluster_ptrs, cluster_data = self.cluster_obj.create_clusters(
                round(resolution, 3)
            )

            allele_data = self.processing_cluster_data(
                cluster_data, cluster_ptrs, position_to_allele
            )
            ref_fasta_file = self.save_reference_alleles(
                allele_data["reference_alleles"]
            )

            # evaluate clusters aginst blast results
            stderr.print(f"Evaluating clusters for {self.locus_name}")
            evaluation_obj = taranis.eval_cluster.EvaluateCluster(
                self.fasta_file, self.locus_name, self.output
            )
            evaluation_result = evaluation_obj.evaluate_clusters(
                allele_data["alleles_in_cluster"],
                allele_data["cluster_data"],
                ref_fasta_file,
            )
            if evaluation_result["result"] == "OK" or resolution >= 1:
                # delete blast database used for evaluation
                _ = evaluation_obj.delete_blast_db_folder()
                break
            stderr.print(
                f"[yellow]{self.locus_name} resolution {resolution} not good enough. Increasing resolution"
            )
            log.info(
                "%s resolution %s not good enough. Increasing resolution",
                self.locus_name,
                resolution,
            )

        return {
            "cluster_data": allele_data["cluster_data"],
            "evaluation": evaluation_result,
        }


def parallel_execution(
    fasta_file: str,
    output: str,
    eval_cluster: bool,
    kmer_size: int,
    sketch_size: int,
    cluster_resolution: float,
    seed: int,
):
    """Parallel execution of the reference alleles creation

    Args:
        fasta_file (str): file name included path for locus
        output (str): output folder
        eval_cluster (bool): True if cluster evaluation must be done
        kmer_size (int): kmer size for mash distance
        sketch_size (int): sketch size for mash distance
        cluster_resolution (float): resolution for clustering
        seed (int): seed for random number generator
    """
    ref_alleles_obj = taranis.reference_alleles.ReferenceAlleles(
        fasta_file,
        output,
        eval_cluster,
        kmer_size,
        sketch_size,
        cluster_resolution,
        seed,
    )
    return ref_alleles_obj.create_ref_alleles()


def collect_statistics(data_alleles: list, eval_cluster: bool, out_folder: str) -> None:
    """Collect the individual statistics for each locus to create graphics

    Args:
        data_alleles (list): list having two dictionaries, cluster_data for
            information and evalluation for the result of evaluating
        eval_cluser (bool): True if evaluation data exists to dump this info
        out_folder (str): folder to save graphics
    """

    def stats_graphics(stats_folder: str, cluster_alleles: dict) -> None:
        """Create the statistics graphics. Bar graphic for number of cluster
            per alleles

        Args:
            stats_folder (str): folder path to store graphic
            cluster_alleles (dict): contain number of clusters as key and number
                of alleles having the same cluster number as value
        """
        stderr.print("Creating graphics")
        log.info("Creating graphics")
        graphic_folder = os.path.join(stats_folder, "graphics")
        _ = taranis.utils.create_new_folder(graphic_folder)
        cluster, alleles = zip(*cluster_alleles.items())
        _ = taranis.utils.create_graphic(
            graphic_folder,
            "num_genes_per_allele.png",
            "bar",
            cluster,
            alleles,
            ["Gene", "Number of clusters"],
            "Number of cluster per gene",
        )

    # split into cluster_data and evaluation_data
    cluster_data = []
    eval_data = []
    cluster_data_graph = {}
    clusters_list = []
    stderr.print("Process starts for collecting statistics")
    log.info("Process starts for collecting statistics")
    # split the data into cluster and evaluation
    for d_allele in data_alleles:
        cluster_data.append(d_allele["cluster_data"])
        eval_data.append(d_allele["evaluation"])
    # collect the number of clusters for each allele
    for c_data in cluster_data:
        cluster_number = len(c_data)
        # get data for graphic
        cluster_data_graph[cluster_number] = (
            cluster_data_graph.get(cluster_number, 0) + 1
        )
        # collect cluster information
        for c_idx, c_value in dict(sorted(c_data.items())).items():
            clusters_list.append(
                c_value["locus_name"]
                + ","
                + str(c_idx)
                + ","
                + str(round(c_value["avg"], 2))
                + ","
                + c_value["center_id"]
                + ","
                + str(c_value["n_seq"])
            )
    heading = "Locus name,cluster number,average,center allele,number of sequences"
    summary_file = os.path.join(out_folder, "evaluate_cluster", "cluster_summary.csv")
    locus_clustering_file = os.path.join(
        out_folder, "evaluate_cluster", "cluster_per_locus.csv"
    )
    with open(summary_file, "w") as fo:
        fo.write(heading + "\n")
        fo.write("\n".join(clusters_list) + "\n")

    _ = stats_graphics(out_folder, cluster_data_graph)

    with open(locus_clustering_file, "w") as fo:
        fo.write("number of clusters, number of locus\n")
        sorted_clusters = dict(sorted(cluster_data_graph.items()))
        for key, value in sorted_clusters.items():
            fo.write(str(key) + "," + str(value) + "\n")

    if eval_cluster:
        heading = "Locus name,cluster number,result,alleles not match in blast,alleles not found in cluster"
        eval_file = os.path.join(
            out_folder, "evaluate_cluster", "cluster_evaluation.csv"
        )
        with open(eval_file, "w") as fo:
            fo.write(heading + "\n")
            for eval in eval_data:
                fo.write("\n".join(eval["individual"]) + "\n")

    return
