import argparse
import glob
import json
import os
from difflib import SequenceMatcher

import mantel
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy


def fill_triangle_matrix(mash_tabpath):
    with open(mash_tabpath, "r") as file:
        mashvals = [list(map(float, line.split())) for i, line in enumerate(file) if i > 0]
    matrix_size = len(mashvals)
    for i in range(matrix_size):
        for j in range(i + 1):
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
        if condensed_mash.shape[0] <= 3:
            print(f"Locus in file {blast_filename} has less than 3 alleles.\n")
            failed_tabs.append(blast_tabpath)
            continue
        try:
            condensed_blast = scipy.spatial.distance.squareform(
                inverted_blast, force="tovector", checks=True
            )
        except ValueError:
            print(f"{blast_tabpath} is not symmetric, skipped")
            failed_tabs.append(blast_tabpath)
            continue
        permutations = int(1 / pval)
        if condensed_mash.shape != condensed_blast.shape:
            print("Blast and mash matrizes have different shape.\n")
            failed_tabs.append(blast_tabpath)
            continue
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
    print(f"{len(failed_tabs)} blast matrixes could not be analyzed due to non-symmetrical, less than three alleles or different shape: {failed_tabs}")
    return mantel_summary


# Argument parser setup
parser = argparse.ArgumentParser(description='Process the Mantel test for genetic data.')
parser.add_argument('root_path', type=str, help='The root directory containing the datasets.')
args = parser.parse_args()

# Use the root_path argument
root_path = args.root_path

datasets = ["brucella", "listeria", "mtuberculosis"]
all_results = {}

for dataset in datasets:
    blast_paths = sorted(glob.glob(os.path.join(root_path, dataset, "blast", "*.csv")))
    mash_paths = sorted(glob.glob(os.path.join(root_path, dataset, "mash", "*.txt")))
    mantel_summary = mantel_tester(blast_paths, mash_paths, pval=0.01)  
    all_results[dataset] = mantel_summary
    with open(f"mantel_test_pval001_{dataset}.json", "w") as f:
        json.dump(mantel_summary, f)

# Create DataFrame for visualization
results_series = pd.Series(
    {
        (dataset, key): value["veridical_correlation"]
        for dataset, results in all_results.items()
        for key, value in results.items()
    }
)

results_df = pd.DataFrame(results_series).reset_index()
results_df.columns = ["Dataset", "Locus", "Veridical Correlation"]

# Crea el boxplot
fig, ax = plt.subplots(figsize=(10, 6))  # Dimensiones en pulgadas (ancho, alto)
results_df.boxplot(by="Dataset", column=["Veridical Correlation"], grid=False, ax=ax)
plt.title("Dataset mash-blast correlation")  # Título opcional
plt.suptitle("")  # Elimina el título por defecto
plt.xlabel("Dataset")  # Etiqueta para el eje x
plt.ylabel("Mantel correlation value")  # Etiqueta para el eje y

# Guarda el boxplot como PNG
plt.savefig(
    "boxplot_mantel_test.png", dpi=300, bbox_inches="tight"
)  # Guarda con alta resolución y ajusta el borde
plt.close()  # Cierra la figura para liberar memoria
