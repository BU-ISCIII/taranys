import numpy as np
import pandas as pd
import glob
import os
import mantel
import scipy
from difflib import SequenceMatcher
import json
import matplotlib.pyplot as plt


def fill_triangle_matrix(mash_tabpath):
    with open(mash_tabpath, "r") as file:
        mashvals = [list(map(float, line.split())) for line in file]

    matrix_size = len(mashvals)

    for i in range(matrix_size):
        for j in range(i + 1):  # Only fill values up to the diagonal
            mashvals[i][j] = mashvals[i][j]
    full_mashtab = mashvals
    tri_mashtable = pd.DataFrame(full_mashtab).fillna(0)
    tri_mashtable_clean = tri_mashtable.drop(tri_mashtable.columns[0], axis=1)
    tri_mashtable_clean[tri_mashtable_clean.columns[-1] + 1] = float(0)
    tri_masharray = tri_mashtable_clean.values
    masharray_transraw = tri_masharray.T
    masharray_clean = np.nan_to_num(masharray_transraw, nan=0.0)
    masharray_full = (
        tri_masharray + masharray_transraw - np.diag(np.diag(masharray_clean))
    )
    return masharray_full


def take_upper_tri_and_dup(full_dist_matrix):
    upper_triangle_matrix = np.triu(full_dist_matrix)
    full_matrix = (
        upper_triangle_matrix
        + upper_triangle_matrix.T
        - np.diag(np.diag(upper_triangle_matrix))
    )
    return full_matrix


def mantel_tester(blast_paths, mash_paths, pval=0.01):
    mantel_summary = {}
    failed_tabs = []
    for blast_tabpath, mash_tabpath in zip(blast_paths, mash_paths):
        blast_filename = os.path.basename(blast_tabpath)
        mash_filename = os.path.basename(mash_tabpath)
        match = SequenceMatcher(
            None, blast_filename, mash_filename
        ).find_longest_match()
        common_name = blast_filename[match.a : match.a + match.size].strip(".")

        blastable = pd.read_csv(blast_tabpath)
        blastarray = blastable.drop(blastable.columns[0], axis=1).to_numpy()
        mirror_blastarray = take_upper_tri_and_dup(blastarray)
        inverted_blast = 100 - mirror_blastarray

        masharray_full = fill_triangle_matrix(mash_tabpath)

        condensed_mash = scipy.spatial.distance.squareform(
            masharray_full, force="tovector", checks=True
        )
        try:
            condensed_blast = scipy.spatial.distance.squareform(
                inverted_blast, force="tovector", checks=True
            )
        except ValueError:
            print(f"{blast_tabpath} is not symmetric, skipped")
            failed_tabs.append(blast_tabpath)
            continue
        permutations = int(1 / pval)
        result = mantel.test(
            condensed_mash, condensed_blast, perms=permutations, method="pearson"
        )
        print(
            f"Results from mantel test between {blast_filename} and {mash_filename}:",
            f"veridical-correlation = {result.r} | p-value = {result.p}",
        )
        mantel_summary[common_name] = {
            "veridical_correlation": result.r,
            "p_value": result.p,
            "z_score": result.z,
        }
    print(f"{len(failed_tabs)} blast matrixes where non-symmetrical: {failed_tabs}")
    return mantel_summary


blast_paths = sorted(glob.glob("blast/*.csv"))
mash_paths = sorted(glob.glob("mash/mash*.txt"))
mantel_summary = mantel_tester(blast_paths, mash_paths, pval=0.01)

with open("mantel_test_pval001.json", "w") as f:
    json.dump(mantel_summary, f)

mantel_df_pval001 = pd.DataFrame.from_dict(mantel_summary)
mantel_df_pval001_tr = mantel_df_pval001.T
mantel_df_pval001_tr.boxplot(column=["veridical_correlation"], return_type="axes")
plt.show()
