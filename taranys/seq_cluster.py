"""
Code extracted from ALPHAClust on 03/28/2024.
Chiu, J.K.H., Ong, R.TH. Clustering biological sequences with dynamic sequence similarity threshold.
BMC Bioinformatics 23, 108 (2022). https://doi.org/10.1186/s12859-022-04643-9
https://github.com/phglab/ALFATClust/blob/main/main/modules/SeqCluster.py
"""

from itertools import combinations, combinations_with_replacement

import igraph as ig
import leidenalg
import numpy as np
from scipy.sparse import coo_matrix


class SeqCluster:
    @classmethod
    def __init__(cls, res_start, seed, is_verbose=True):
        cls._res_param_start = res_start
        cls._res_param_end = 0.99
        cls._res_param_step_size = 0.025
        cls._precision = 3
        cls._seed = seed
        cls._is_verbose = is_verbose

    @classmethod
    def disable_verbose(cls):
        cls._is_verbose = False

    @staticmethod
    def _convert_to_index_pos(mtrx_idxs):
        idx_type = mtrx_idxs.dtype

        if idx_type == bool:
            return np.argwhere(mtrx_idxs).flatten()
        elif idx_type == np.int:
            return mtrx_idxs
        else:
            return None

    @classmethod
    def _cal_cluster_avg_edge_weight(
        cls,
        global_edge_weight_mtrx,
        src_cluster_edge_counts,
        src_cluster_row_idxs,
        src_cluster_col_idxs=None,
    ):
        if src_cluster_col_idxs is None:
            src_cluster_col_idxs = src_cluster_row_idxs

        src_cluster_mtrx_idxs = np.ix_(src_cluster_row_idxs, src_cluster_col_idxs)

        if np.any(global_edge_weight_mtrx[src_cluster_mtrx_idxs] < 0):
            return -1

        if src_cluster_edge_counts is None:
            row_idx_pos = cls._convert_to_index_pos(src_cluster_row_idxs)
            col_idx_pos = cls._convert_to_index_pos(src_cluster_col_idxs)
            num_of_diag_elements = np.intersect1d(row_idx_pos, col_idx_pos).size
            num_of_non_diag_elements = (
                row_idx_pos.size * col_idx_pos.size - num_of_diag_elements
            )
            if num_of_non_diag_elements == 0:
                return 1

            return (
                np.sum(global_edge_weight_mtrx[src_cluster_mtrx_idxs])
                - num_of_diag_elements
            ) / num_of_non_diag_elements
        else:
            return np.average(
                global_edge_weight_mtrx[src_cluster_mtrx_idxs],
                weights=np.tril(src_cluster_edge_counts[src_cluster_mtrx_idxs]),
            )

    @classmethod
    def _update_edge_weight_mtrx(
        cls,
        src_cluster_ptrs,
        global_edge_weight_mtrx,
        avg_intra_super_cluster_edge_weights,
        super_cluster_pairs_to_isolate,
        src_cluster_edge_counts=None,
    ):
        row_idxs = list()
        col_idxs = list()
        new_edge_weights = list()
        num_of_super_clusters = np.max(src_cluster_ptrs) + 1

        for super_cluster_pair in combinations_with_replacement(
            range(num_of_super_clusters), 2
        ):
            if super_cluster_pair in super_cluster_pairs_to_isolate:
                continue

            if super_cluster_pair[0] == super_cluster_pair[1]:
                row_idxs.append(super_cluster_pair[0])
                col_idxs.append(super_cluster_pair[0])
                new_edge_weights.append(
                    avg_intra_super_cluster_edge_weights[super_cluster_pair[0]]
                )
                continue

            src_cluster_bool_ptrs1 = src_cluster_ptrs == super_cluster_pair[0]
            src_cluster_bool_ptrs2 = src_cluster_ptrs == super_cluster_pair[1]
            inter_super_cluster_edge_weight = cls._cal_cluster_avg_edge_weight(
                global_edge_weight_mtrx,
                src_cluster_edge_counts,
                src_cluster_bool_ptrs1,
                src_cluster_bool_ptrs2,
            )

            if inter_super_cluster_edge_weight <= 0:
                continue

            row_idxs += super_cluster_pair
            col_idxs += [super_cluster_pair[1], super_cluster_pair[0]]
            new_edge_weights += [inter_super_cluster_edge_weight] * 2

        global_edge_weight_mtrx = coo_matrix(
            (np.array(new_edge_weights), (np.array(row_idxs), np.array(col_idxs))),
            shape=(num_of_super_clusters, num_of_super_clusters),
        ).toarray()

        global_edge_weight_mtrx[global_edge_weight_mtrx == 0] = (
            -1 * global_edge_weight_mtrx.size
        )

        return global_edge_weight_mtrx

    @staticmethod
    def _convert_cluster_ptrs(src_cluster_ptrs, last_seq_cluster_ptrs):
        if src_cluster_ptrs.size == last_seq_cluster_ptrs.size:
            return src_cluster_ptrs

        output_seq_cluster_ptrs = np.array([-1] * last_seq_cluster_ptrs.size)

        for cluster_id in range(np.max(src_cluster_ptrs) + 1):
            src_cluster_ids = np.argwhere(src_cluster_ptrs == cluster_id).flatten()
            output_seq_cluster_ptrs[np.isin(last_seq_cluster_ptrs, src_cluster_ids)] = (
                cluster_id
            )

        return output_seq_cluster_ptrs

    @staticmethod
    def _count_intra_cluster_edges(seq_cluster_ptrs):
        row_idxs = list()
        col_idxs = list()
        intra_cluster_edge_counts = list()
        num_of_clusters = np.max(seq_cluster_ptrs) + 1

        for cluster_pair in combinations_with_replacement(range(num_of_clusters), 2):
            cluster_size1 = np.count_nonzero(seq_cluster_ptrs == cluster_pair[0])
            if cluster_pair[0] == cluster_pair[1]:
                row_idxs.append(cluster_pair[0])
                col_idxs.append(cluster_pair[0])
                intra_cluster_edge_counts.append(
                    int(cluster_size1 * (cluster_size1 - 1) / 2)
                )
            else:
                row_idxs += cluster_pair
                col_idxs += [cluster_pair[1], cluster_pair[0]]
                edge_count = cluster_size1 * np.count_nonzero(
                    seq_cluster_ptrs == cluster_pair[1]
                )
                intra_cluster_edge_counts += [edge_count, edge_count]

        return coo_matrix(
            (
                np.array(intra_cluster_edge_counts),
                (np.array(row_idxs), np.array(col_idxs)),
            ),
            shape=(num_of_clusters, num_of_clusters),
        ).toarray()

    @staticmethod
    def _sort_src_clusters(
        global_edge_weight_mtrx, src_cluster_ids, src_cluster_edge_counts
    ):
        cluster_mtrx_idxs = np.ix_(src_cluster_ids, src_cluster_ids)
        all_sorted_edge_row_idxs, all_sorted_edge_col_idxs = np.unravel_index(
            np.argsort(global_edge_weight_mtrx[cluster_mtrx_idxs], axis=None),
            shape=(src_cluster_ids.size, src_cluster_ids.size),
        )
        all_sorted_edge_idxs = np.array(
            list(zip(all_sorted_edge_row_idxs, all_sorted_edge_col_idxs))
        )

        sorted_src_cluster_ids = list()
        proc_src_cluster_ids = set()

        for src_cluster_idx1, src_cluster_idx2 in np.flipud(
            all_sorted_edge_idxs[
                all_sorted_edge_idxs[:, 0] < all_sorted_edge_idxs[:, 1]
            ]
        ):
            src_cluster_id1 = src_cluster_ids[src_cluster_idx1]
            src_cluster_id2 = src_cluster_ids[src_cluster_idx2]

            if src_cluster_id1 in proc_src_cluster_ids:
                if src_cluster_id2 in proc_src_cluster_ids:
                    continue

                sorted_src_cluster_ids.append(src_cluster_id2)
                proc_src_cluster_ids.add(src_cluster_id2)
                continue
            elif src_cluster_id2 in proc_src_cluster_ids:
                sorted_src_cluster_ids.append(src_cluster_id1)
                proc_src_cluster_ids.add(src_cluster_id1)
                continue

            if (
                src_cluster_edge_counts[src_cluster_id1, src_cluster_id1]
                < src_cluster_edge_counts[src_cluster_id2, src_cluster_id2]
            ):
                sorted_src_cluster_ids.append(src_cluster_id2)
                proc_src_cluster_ids.add(src_cluster_id2)
                sorted_src_cluster_ids.append(src_cluster_id1)
                proc_src_cluster_ids.add(src_cluster_id1)
            else:
                sorted_src_cluster_ids.append(src_cluster_id1)
                proc_src_cluster_ids.add(src_cluster_id1)
                sorted_src_cluster_ids.append(src_cluster_id2)
                proc_src_cluster_ids.add(src_cluster_id2)

        return sorted_src_cluster_ids

    @classmethod
    def _bin_src_clusters(
        cls, src_cluster_ids, global_edge_weight_mtrx, src_cluster_edge_counts
    ):
        src_cluster_bins = list()
        avg_intra_bin_edge_weights = dict()

        for src_cluster_id in cls._sort_src_clusters(
            global_edge_weight_mtrx, src_cluster_ids, src_cluster_edge_counts
        ):
            if len(src_cluster_bins) == 0:
                src_cluster_bins.append([src_cluster_id])
                avg_intra_bin_edge_weights[str(src_cluster_id)] = (
                    global_edge_weight_mtrx[src_cluster_id, src_cluster_id]
                )
                continue

            best_bin = None
            max_merge_score = 0

            for cluster_bin in src_cluster_bins:
                bin_to_ext_src_cluster_avg_edge_weight = (
                    cls._cal_cluster_avg_edge_weight(
                        global_edge_weight_mtrx,
                        src_cluster_edge_counts,
                        cluster_bin,
                        [src_cluster_id],
                    )
                )
                merge_score = (
                    2 * bin_to_ext_src_cluster_avg_edge_weight
                    - cls._res_param_end
                    - avg_intra_bin_edge_weights["-".join(map(str, cluster_bin))]
                )
                if merge_score > max_merge_score:
                    best_bin = cluster_bin
                    max_merge_score = merge_score

            if best_bin is None:
                src_cluster_bins.append([src_cluster_id])
                avg_intra_bin_edge_weights[str(src_cluster_id)] = (
                    global_edge_weight_mtrx[src_cluster_id, src_cluster_id]
                )
            else:
                del avg_intra_bin_edge_weights["-".join(map(str, best_bin))]
                best_bin.append(src_cluster_id)
                avg_intra_bin_edge_weights["-".join(map(str, best_bin))] = (
                    cls._cal_cluster_avg_edge_weight(
                        global_edge_weight_mtrx, src_cluster_edge_counts, best_bin
                    )
                )

        return list(map(np.array, src_cluster_bins)), avg_intra_bin_edge_weights

    @classmethod
    def _verify_clusters(
        cls,
        candidate_src_cluster_ptrs,
        global_edge_weight_mtrx,
        src_cluster_edge_counts,
        is_first_iter,
    ):
        super_cluster_pairs_to_isolate = set()
        avg_intra_super_cluster_edge_weights = dict()

        if is_first_iter:
            src_cluster_ptrs = candidate_src_cluster_ptrs

            for cluster_id in range(np.max(src_cluster_ptrs) + 1):
                src_cluster_bool_ptrs = src_cluster_ptrs == cluster_id
                avg_intra_super_cluster_edge_weights[cluster_id] = (
                    cls._cal_cluster_avg_edge_weight(
                        global_edge_weight_mtrx,
                        src_cluster_edge_counts,
                        src_cluster_bool_ptrs,
                    )
                )
        else:
            src_cluster_ptrs = np.full(candidate_src_cluster_ptrs.size, -1)

            for cluster_id in range(np.max(candidate_src_cluster_ptrs) + 1):
                src_cluster_ids = np.argwhere(
                    candidate_src_cluster_ptrs == cluster_id
                ).flatten()

                if src_cluster_ids.size == 1:
                    assigned_super_cluster_id = np.max(src_cluster_ptrs) + 1
                    src_cluster_ptrs[src_cluster_ids] = assigned_super_cluster_id
                    avg_intra_super_cluster_edge_weights[assigned_super_cluster_id] = (
                        global_edge_weight_mtrx[src_cluster_ids[0], src_cluster_ids[0]]
                    )
                    continue

                (
                    qual_src_cluster_bins,
                    avg_intra_bin_edge_weights,
                ) = cls._bin_src_clusters(
                    src_cluster_ids,
                    global_edge_weight_mtrx,
                    src_cluster_edge_counts,
                )

                assigned_super_cluster_ids = list()
                for src_cluster_bin in qual_src_cluster_bins:
                    assigned_super_cluster_ids.append(np.max(src_cluster_ptrs) + 1)
                    src_cluster_ptrs[src_cluster_bin] = assigned_super_cluster_ids[-1]
                    avg_intra_super_cluster_edge_weights[
                        assigned_super_cluster_ids[-1]
                    ] = avg_intra_bin_edge_weights["-".join(map(str, src_cluster_bin))]

                if len(assigned_super_cluster_ids) > 1:
                    super_cluster_pairs_to_isolate |= set(
                        combinations(assigned_super_cluster_ids, 2)
                    )

        return (
            src_cluster_ptrs,
            avg_intra_super_cluster_edge_weights,
            super_cluster_pairs_to_isolate,
        )

    @classmethod
    def cluster_seqs(cls, global_edge_weight_mtrx):
        num_of_seqs = global_edge_weight_mtrx.shape[0]

        last_seq_cluster_ptrs = np.arange(num_of_seqs)

        global_edge_weight_mtrx[global_edge_weight_mtrx == 0] = (
            -1 * global_edge_weight_mtrx.size
        )

        comm_graph = ig.Graph.Weighted_Adjacency(
            global_edge_weight_mtrx.tolist(), mode=1, loops=False
        )
        res_param_end = round(
            cls._res_param_end + cls._res_param_step_size, cls._precision
        )

        for res_param in np.arange(
            cls._res_param_start, res_param_end, cls._res_param_step_size
        ):
            res_param = round(res_param, cls._precision)

            graph_partitions = leidenalg.find_partition(
                comm_graph,
                leidenalg.CPMVertexPartition,
                weights="weight",
                n_iterations=-1,
                resolution_parameter=res_param,
                seed=cls._seed,
            )

            candidate_src_cluster_ptrs = np.array(graph_partitions.membership)
            comm_graph = None

            is_first_iter = res_param == cls._res_param_start

            if is_first_iter:
                src_cluster_edge_counts = None
            else:
                src_cluster_edge_counts = cls._count_intra_cluster_edges(
                    last_seq_cluster_ptrs
                )

            (
                src_cluster_ptrs,
                avg_intra_super_cluster_edge_weights,
                super_cluster_pairs_to_isolate,
            ) = cls._verify_clusters(
                candidate_src_cluster_ptrs,
                global_edge_weight_mtrx,
                src_cluster_edge_counts,
                is_first_iter,
            )

            global_edge_weight_mtrx = cls._update_edge_weight_mtrx(
                src_cluster_ptrs,
                global_edge_weight_mtrx,
                avg_intra_super_cluster_edge_weights,
                super_cluster_pairs_to_isolate,
                src_cluster_edge_counts,
            )

            last_seq_cluster_ptrs = cls._convert_cluster_ptrs(
                src_cluster_ptrs, last_seq_cluster_ptrs
            )

            if cls._is_verbose:
                proc_msg = "{} clusters obtained at average estimated similarity {}"
                print(proc_msg.format(np.max(last_seq_cluster_ptrs) + 1, res_param))

            if np.max(last_seq_cluster_ptrs) == 0 or np.all(
                np.triu(global_edge_weight_mtrx, k=1) <= 0
            ):
                if cls._is_verbose:
                    print("No more cluster available for further merging")

                break

            comm_graph = ig.Graph.Weighted_Adjacency(
                global_edge_weight_mtrx.tolist(), mode=1, loops=False
            )

        return last_seq_cluster_ptrs
