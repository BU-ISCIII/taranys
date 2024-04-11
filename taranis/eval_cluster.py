import io
import logging
import numpy as np
import rich.console
import os
import taranis.utils
import taranis.blast
from Bio import SeqIO

log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)


class EvaluateCluster:
    def __init__(self, locus_path: str, locus_name: str, eval_id: float, output: str):
        """EvaluateCluster instance creation

        Args:
            locus_path (str): path of the locus
            locus_name (str): locus name
            output (str): folder to store results
        """
        self.locus_path = locus_path
        self.locus_name = locus_name
        self.eval_id = eval_id

        self.output = os.path.join(output, "evaluate_cluster")
        taranis.utils.create_new_folder(self.output)
        # locus_blast_dir = os.path.join(self.output, locus_name)
        self.blast_obj = taranis.blast.Blast("nucl")
        _ = self.blast_obj.create_blastdb(locus_path, self.output)
        return

    def delete_blast_db_folder(self):
        """Delete blast db folder"""
        taranis.utils.delete_folder(os.path.join(self.output, self.locus_name))

    def find_cluster_from_ref_allele(self, cluster_ref_alleles: dict) -> dict:
        """Create a dictionary to map de cluster belongs to the reference allele

        Args:
            cluster_ref_alleles (dict): values collected for statistics with the
                cluster id and the reference allele name

        Returns:
            dict: relation between reference allele name and the cluster
        """
        return dict(
            [(value["center_id"], c_id) for c_id, value in cluster_ref_alleles.items()]
        )

    def summary(self, cluster_data: dict) -> dict:
        """Create the summary information from the individual result for each
           cluster

        Args:
            cluster_data (dict): cluster evaluation

        Returns:
            dict: summary table for getting nice presentation of evaluation data
                and global result for the locus
        """
        summary_table = []
        summary_data = {"result": "OK"}
        # heading = "Locus name,result,alleles not found,alleles not in cluster"
        # summary_table.append(heading)
        sorted_cluster = sorted(cluster_data.keys())
        for cluster_id in sorted_cluster:
            row_data = [self.locus_name, str(cluster_id)]
            row_data.append(cluster_data[cluster_id]["result"])
            row_data.append(
                ";".join(cluster_data[cluster_id]["alleles_not_found"])
                if "alleles_not_found" in cluster_data[cluster_id]
                else "-"
            )
            row_data.append(
                ";".join(cluster_data[cluster_id]["alleles_not_in_cluster"])
                if "alleles_not_in_culster" in cluster_data[cluster_id]
                else "-"
            )
            if cluster_data[cluster_id]["result"] == "NOK":
                summary_data["result"] = "NOK"
            summary_table.append(",".join(row_data))
        summary_data["individual"] = summary_table
        return summary_data

    def validate_cluster(self, blast_result: list, cluster_data: list) -> dict:
        """For cluster validation, the sequence id matched in blast are compared
            with the cluster sequences. Return False validation if there are
            difference between them.

        Args:
            blast_result (list): blast matches results
            cluster_data (list): allele names for the cluster to evaluate

        Returns:
            dict: result of the evaluation
        """
        # index of sequence id
        sseqid = 1
        blast_alleles = []
        alleles_not_in_cluster = []
        for match in blast_result:
            blast_allele = match.split("\t")[sseqid]
            if blast_allele in cluster_data:
                blast_alleles.append(blast_allele)
            else:
                alleles_not_in_cluster.append(blast_allele)

        if len(cluster_data) == len(set(blast_alleles)):
            return {"validation": True}
        result = {"validation": False}
        # convert list to numpy array to find out differences
        c_alleles_np = np.array(list(cluster_data))
        blast_alleles_np = np.array(blast_alleles)
        result["alleles_not_found"] = np.setdiff1d(
            c_alleles_np, blast_alleles_np
        ).tolist()
        result["alleles_not_in_cluster"] = np.setdiff1d(
            blast_alleles_np, c_alleles_np
        ).tolist()
        return result

    def evaluate_clusters(
        self, cluster_alleles: dict, cluster_ref_alleles: dict, ref_alleles_file: str
    ) -> list:
        """Perform clusted evaluation comparing for each clusted defined in
            previous step with searching the matches that blast found running
            witha 90% of percentage of identity

        Args:
            cluster_alleles (dict): contains the cluster id as dict and the list
                of allele names as value
            cluster_ref_alleles (dict): statistics information for each cluster
                to fetch the reference allele for each cluster
            ref_alleles_file (str): reference alleles to get the seqence for the
                reference allele

        Returns:
            list: evaluation imformation for each cluster
        """
        reference_alleles = {}
        evaluation_alleles = {}
        ref_allele_in_cluster = self.find_cluster_from_ref_allele(cluster_ref_alleles)
        with open(ref_alleles_file, "r") as fh:
            for record in SeqIO.parse(fh, "fasta"):
                reference_alleles[record.id] = str(record.seq)

        for r_id, r_seq in reference_alleles.items():
            # create file in memory to increase speed
            query_file = io.StringIO()
            query_file.write(">" + r_id + "\n" + r_seq)
            query_file.seek(0)
            blast_result = self.blast_obj.run_blast(
                query_file.read(), perc_identity=self.eval_id, query_type="stdin"
            )
            # Close object and discard memory buffer
            query_file.close()

            cluster_id = ref_allele_in_cluster[r_id]
            result_eval = self.validate_cluster(
                blast_result, cluster_alleles[cluster_id]
            )
            evaluation_alleles[cluster_id] = {}
            if result_eval["validation"] is False:
                evaluation_alleles[cluster_id]["result"] = "NOK"
                if len(result_eval["alleles_not_found"]) > 0:
                    evaluation_alleles[cluster_id]["alleles_not_found"] = result_eval[
                        "alleles_not_found"
                    ]
                if len(result_eval["alleles_not_in_cluster"]) > 0:
                    evaluation_alleles[cluster_id]["alleles_not_in_cluster"] = (
                        result_eval["alleles_not_in_cluster"]
                    )
            else:
                evaluation_alleles[cluster_id]["result"] = "OK"
        return self.summary(evaluation_alleles)
