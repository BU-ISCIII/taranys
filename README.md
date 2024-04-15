# Taranis

- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Installation](#installation)
  - [Install from source](#install-from-source)
  - [Install using conda](#install-using-conda)
- [Quick usage](#quick-usage)
- [Usage](#usage)
- [Output](#output)
- [Illustrated pipeline](#illustrated-pipeline)

## Introduction

**Taranis** is a computational stand-alone pipeline for **gene-by-gene allele calling analysis** based on BLASTn using  whole genome (wg) and core genome (cg) multilocus sequence typing (MLST) schemas on complete or draft genomes resulting from de novo assemblers, while tracking helpful and informative data among the process.

Taranis includes four main functionalities: MLST **schema analysis**, gene-by-gene **allele calling**, **reference alleles** obtainment for allele calling analysis and the final **distance matrix** construction.

## Dependencies

- Python >=3.8
- NCBI_blast >= v2.9
- prokka >=1.14.6
- mafft = 7.520
- mash >=2
- biopython v1.81
- pandas v2.1.1
- plotly v5.17.0
- numpy v1.26.0

## Installation

#### Install from source

Install all dependencies and add them to $PATH.

`git clone https://github.com/BU-ISCIII/taranis.git`

Add taranis and ./bin to $PATH.

#### Install using conda

This option is recomended.

Install Anaconda3.

`conda install -c conda-forge -c bioconda -c defaults taranis`

Wait for the environment to solve. <br>
Ignore warnings/errors.

## Quick usage

- **analyze_schema mode:**

  Schema analysis:

```
taranis analyze_schema \
-inputdir schema_dir \
-output output_analyze_schema_dir
--ouput-allele-annotation annotation_dir
```

  Schema analysis for removing duplicated, subsequences and no CDS alleles:

```
taranis analyze_schema \
-inputdir schema_dir \
-output output_analyze_schema_dir \
--remove-subsets \
--remove-duplicated \
--remove-no-cds \
--ouput-allele-annotation annotation_dir \
--genus prokka_genus_name \
--usegenus prokka genus-specific BLAST database \
--species prokka_species_name \
--cpus number_of_cpus

```

- **reference_alleles mode:**

  Get reference alleles:

```
taranis reference_alleles \
-s schema_dir \
-o output_reference_alleles_dir \
--eval-cluster \
--cpus number_of_cpus \
--force overwrite output dir
```

  Reference alleles with clustering settings:

```
taranis reference_alleles \
-s schema_dir \
-o output_reference_alleles_dir \
--eval-cluster \
-k k-mer size for mash \
-S Sketch size for mash \
-r resolution used for clustering \
--cpus number_of_cpus \
--force overwrite output dir
```

- **allele_calling mode:**

  Run allele calling:

```
taranis allele_calling \
-s schema_dir \
-a annotation_file \
-r reference_alleles_dir \
-o output_allele_calling_dir \
-t threshold to consider in blast \
-p percentage of identity to consider in blast \
-q threshold to consider as TPR \
-i increase number of nucleotides to find stop codon \
--snp Create SNP file \
--cpus number_of_cpus \
--alignment Create aligment files \
samples_dir 
```

  Allele calling for blast and threshold settings:

```
taranis allele_calling \
-s schema_dir \
-a annotation_file \
-r reference_alleles_dir \
-o output_allele_calling_dir \
-t threshold to consider in blast \
-p percentage of identity to consider in blast \
-q threshold to consider as TPR \
-i increase number of nucleotides to find stop codon \
--snp Create SNP file \
--cpus number_of_cpus \
--alignment Create aligment files \
samples_dir 
```

- **distance_matrix mode:**

  Get distance matrix:

```
taranis distance_matrix \
-a allele_calling_match.csv file \
-o distance_matrix_dir
--force overwrite output folder
```

Distance matrix with threshold settings:

```
taranis distance_matrix \
-a allele_calling_match.csv file \
-o distance_matrix_dir
-l threshold for missing locus \
-s threshold for missing samples \
--paralog-filter \
--lnf-filter \
--plot-filter \
--force overwrite output folder
```

## Usage

- **analyze_schema mode:**

```
Usage: taranis analyze-schema [OPTIONS]

Options:
  -i, --inputdir PATH             Directory where the schema with the core
                                  gene files are located.   [required]
  -o, --output PATH               Output folder to save analyze schema
                                  [required]
  --remove-subset / --no-remove-subset
                                  Remove allele subsequences from the schema.
  --remove-duplicated / --no-remove-duplicated
                                  Remove duplicated subsequences from the
                                  schema.
  --remove-no-cds / --no-remove-no-cds
                                  Remove no CDS alleles from the schema.
  --output-allele-annot / --no-output-allele-annot
                                  output prokka/allele annotation for all
                                  alleles in locus
  --genus TEXT                    Genus name for Prokka schema genes
                                  annotation. Default is Genus.
  --species TEXT                  Species name for Prokka schema genes
                                  annotation. Default is species
  --usegenus TEXT                 Use genus-specific BLAST databases for
                                  Prokka schema genes annotation (needs
                                  --genus). Default is False.
  --cpus INTEGER                  Number of cpus used for execution
  --help                          Show this message and exit.
```

- **reference_alleles mode:**

```
Usage: taranis reference-alleles [OPTIONS]

Options:
  -s, --schema PATH               Directory where the schema with the core
                                  gene files are located.   [required]
  -o, --output PATH               Output folder to save reference alleles
                                  [required]
  --eval-cluster / --no-eval-cluster
                                  Evaluate if the reference alleles match
                                  against blast with a 90% identity
  -k, --kmer-size INTEGER         Mash parameter for K-mer size.
  -S, --sketch-size INTEGER       Mash parameter for Sketch size
  -r, --cluster-resolution FLOAT  Resolution value used for clustering.
  --seed INTEGER                  Seed value for clustering
  --cpus INTEGER                  Number of cpus used for execution
  --force / --no-force            Overwrite the output folder if it exists
  --help                          Show this message and exit.
```

- **allele_calling mode:**

```
Usage: taranis allele-calling [OPTIONS] ASSEMBLIES...

Options:
  -s, --schema PATH               Directory where the schema with the core
                                  gene files are located.   [required]
  -r, --reference PATH            Directory where the schema reference allele
                                  files are located.   [required]
  -a, --annotation PATH           Annotation file.   [required]
  -t, --threshold FLOAT           Threshold value to consider in blast. Values
                                  from 0 to 1. default 0.8
  -p, --perc-identity INTEGER     Percentage of identity to consider in blast.
                                  default 90
  -o, --output PATH               Output folder to save reference alleles
                                  [required]
  --force / --no-force            Overwrite the output folder if it exists
  --snp / --no-snp                Create SNP file for alleles in assembly in
                                  relation with reference allele
  --alignment / --no-alignment    Create alignment files
  -q, --proteine-threshold INTEGER
                                  Threshold of protein coverage to consider as
                                  TPR. default 90
  -i, --increase-sequence INTEGER
                                  Increase the number of triplet sequences to
                                  find the stop codon. default 20
  --cpus INTEGER                  Number of cpus used for execution
  --help                          Show this message and exit.
```

- **distance_matrix mode:**

```
Usage: taranis distance-matrix [OPTIONS]

Options:
  -a, --alleles PATH              Alleles matrix file from which to obtain
                                  distances between samples  [required]
  -o, --output PATH               Output folder to save distance matrix
                                  [required]
  --force / --no-force            Overwrite the output folder if it exists
  -l, --locus-missing-threshold INTEGER
                                  Threshold for missing alleles in locus,
                                  which loci is excluded from distance matrix
  -s, --sample-missing-threshold INTEGER
                                  Threshold for missing samples, which sample
                                  is excluded from distance matrix
  --paralog-filter / --no-paralog-filter
                                  Consider paralog tags (NIPH, NIPHEM) as
                                  missing values. Default is True
  --lnf-filter / --no-lnf-filter  Consider LNF as missing values. Default is
                                  True
  --plot-filter / --no-plot-filter
                                  Consider PLOT as missing values. Default is
                                  True
  --help                          Show this message and exit.
```

## Output

- **analyze_schema mode:**

  - **FOLDERS and FILES structure:**
  
    - **new_schema**  Contains the new schema.
    - **prokka** Contains the prokka results
    - **statistics** Statistics data
      - **graphics** Plot graphics folder
      - **statistics.csv** Quality statistics showing the following data:

        - allele_name,
        - min_length,
        - max_length,
        - num_alleles,
        - mean_length,
        - good_percent,
        - not a start codon,
        - not a stop codon,
        - Extra in frame stop codon,
        - is not a multiple of three,
        - Duplicate allele,
        - Sub set allele

    - **allele_annotation.tar.gz** Annotation schema file

- **reference_alleles mode:**

  - **FOLDERS and FILES structure:**
  
    - **Clusters** Contains the cluster allele files
      - **[cluster_alleles].txt** cluster allele file
    - **evaluate_cluster**
      - **cluster_evaluuation.csv** Evaluation result with the following info:
        - Locus name
        - cluster number
        - result
        - alleles not match in blast
        - alleles not found in cluster

      - **cluster_per_locus.csv** Number of cluster per locus
        - number of clusters
        - number of locus

      - **cluster_summary.csv** summary data with the following info:
        - Locus name
        - cluster number
        - average
        - center allele
        - number of sequences

    - **graphics** Plot graphics folder
      - **num_genes_per_allele.png** Bar graphic to show the number of clusters per gene

    - **[ref_alleles_locusX].fasta:** One fasta file for each schema locus containing reference alleles for that locus

- **allele_calling mode:**

  - **FOLDERS and FILES structure:**
    - **alignments:** Nucleotide alignment between sequence found in the sample and allele
      - **[locus_name].txt** One file per locus
      - **[locus_name]_multiple_alignment.aln** One file per locus
    - **graphics** Graphics per type of allele classification
      - **ALM_graphic.pnd** Number of ALM in samples.
      - **ASM_graphic.pnd** Number of ASM in samples.
      - **EXEC_graphic.pnd** Number of EXEC in samples.
      - **INF_graphic.pnd** Number of INF in samples.
      - **LNF_graphic.pnd** Number of LNF in samples.
      - **NIPHEM_graphic.pnd** Number of NIPHEM in samples.
      - **NIPH_graphic.pnd** Number of NIPH in samples.
      - **PLOT_graphic.pnd** Number of PLOT in samples.
      - **TPR_graphic.pnd** Number of TPR in samples.
    - **[locus_name]_snp_data** One file per sample
    - **allele_calling_match.csv** Contains the classification for each locus and for all samples
    - **allele_calling_summary.csv** Contains the number of each classification per samples
    - **matching_contig.csv** Summary for each locus in sample with the following data:
      - sample
      - contig
      - core gene
      - reference allele name
      - codification
      - query length
      - match length
      - contig length
      - contig start
      - contig stop
      - direction
      - gene notation
      - product notation
      - reference allele quality
      - protein conversion result
      - match sequence reference
      - allele sequence
      - predicted protein sequence

- **distance_matrix mode:**

  - **FILES:**
    - **filtered_result.tsv:** Filtered allele calling matrix filtered
    - **matrix_distance.tsv:** Samples matrix distance
    - **matrix_distance_filter_report.tsv:** Allele calling matrix filtering report

## Illustrated pipeline
