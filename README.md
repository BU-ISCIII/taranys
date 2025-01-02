# taranys

- [taranys](#taranys)
  - [Introduction](#introduction)
  - [Dependencies](#dependencies)
  - [Installation](#installation)
    - [Install from source](#install-from-source)
    - [Install using conda](#install-using-conda)
  - [Usage](#usage)
    - [**analyze\_schema:**](#analyze_schema)
    - [**reference\_alleles:**](#reference_alleles)
    - [**allele\_calling:**](#allele_calling)
    - [**distance\_matrix:**](#distance_matrix)
  - [Output](#output)

## Introduction

**taranys** is a computational stand-alone pipeline for **gene-by-gene allele calling analysis** based on BLASTn using  whole genome (wg) and core genome (cg) multilocus sequence typing (MLST) schemas on complete or draft genomes resulting from de novo assemblers, while tracking helpful and informative data among the process.

taranys includes four main functionalities: MLST **schema analysis**, gene-by-gene **allele calling**, **reference alleles** obtainment for allele calling analysis and the final **distance matrix** construction.

![taranis_schema](assets/taranis_schema.png)

## Dependencies

- Python >=3.10
- NCBI_blast >= v2.9
- prokka >=1.14.6
- mafft >= 7.505
- mash >= 2
- python deps:
  - igraph>=0.9.8
  - rich>=13.4.1
  - click>=8.1.3
  - leidenalg>=0.9.1
  - questionary>=1.10.0
  - bio>=1.6.0
  - scikit-learn>=1.2.0
  - plotly>=5.11.0
  - kaleido>=0.2.1
  - six>=1.16.0

## Installation

### Install from source

```bash
git clone https://github.com/BU-ISCIII/taranys.git
```

Install dependencies and taranys using conda or micromamba:

```bash
cd /path/to/clonedrepo
micromamba install -f environment.yml
```

### Install using conda

This option is the recommended option for installing taranys.

```bash
micromamba install -c conda-forge -c bioconda -c defaults taranys
```

## Usage

### **analyze_schema:**

To assess the quality of the schema, the following analysis is performed for each locus in the schema:

- The existence of potential duplicate alleles within the same locus is examined.
- The presence of allelic sequences that are partial sequences or subsequences of other alleles is checked.
- Each allele is evaluated to verify whether it is a coding region (CDS) using the translate function from Biopython's Seq class.

The following quality categories are defined:

- “Good” quality: The sequence meets the criteria to be considered a hypothetical CDS.
- “Bad” quality: The sequence fails to meet the criteria for a hypothetical CDS due to one of the following reasons:
  - Lack of a start codon (Bad quality: no start codon).
  - Lack of a stop codon (Bad quality: no stop codon).
  - Simultaneous absence of both start and stop codons (Bad quality: no start stop).
  - The sequence length is not a multiple of 3 (Bad quality: no multiple of three).
  - Presence of multiple stop codons (Bad quality: multiple stop).

Usage:

```bash
Usage: taranys analyze-schema [OPTIONS]

Options:
  -i, --input PATH                Directory where the schema with the core
                                  gene files are located.   [required]
  -o, --output PATH               Output folder to save analyze schema
                                  [required]
  --remove-subset / --no-remove-subset
                                  Remove allele subsequences from the schema.
                                  [default: no-remove-subset]
  --remove-duplicated / --no-remove-duplicated
                                  Remove duplicated subsequences from the
                                  schema.  [default: no-remove-duplicated]
  --remove-no-cds / --no-remove-no-cds
                                  Remove no CDS alleles from the schema.
                                  [default: no-remove-no-cds]
  --output-allele-annot / --no-output-allele-annot
                                  output prokka/allele annotation for all
                                  alleles in locus. Default is True.
                                  [default: output-allele-annot]
  --genus TEXT                    Genus name for Prokka schema genes
                                  annotation. Default is Genus.  [default:
                                  Genus]
  --species TEXT                  Species name for Prokka schema genes
                                  annotation. Default is species  [default:
                                  species]
  --usegenus TEXT                 Use genus-specific BLAST databases for
                                  Prokka schema genes annotation (needs
                                  --genus). Default is False.  [default:
                                  Genus]
  --cpus INTEGER                  Number of cpus used for execution. Default
                                  is 1  [default: 1]
  --help                          Show this message and exit.
```

Example when removing bad quality alleles is not wanted, just statistics outputted:

```bash
taranys analyze-schema \
-input schema_dir \
-output output_analyze_schema_dir
--ouput-allele-annotation annotation_dir
```

Example for removing duplicated, subsequences and no CDS alleles:

```bash
taranys analyze-schema \
-input schema_dir \
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

### **reference_alleles:**

As a preliminary step to the typing analysis, the representative allele or alleles for each locus in the schema are determined. The selected allele will be the one that shows the least dissimilarity compared to all other known alleles for that locus and can be used to detect, with a certain degree of similarity, all alleles in the schema. Leiden algorithm is used for clusting similar sequences.

Usage:

```bash
Usage: taranys reference-alleles [OPTIONS]

Options:
  -s, --schema PATH               Directory where the schema with the core
                                  gene files are located.   [required]
  -o, --output PATH               Output folder to save reference alleles
                                  [required]
  --eval-cluster / --no-eval-cluster
                                  Evaluate if the reference alleles match
                                  against blast with the identity set in eval-
                                  identity param  [default: eval-cluster]
  -k, --kmer-size INTEGER         Mash parameter for K-mer size.  [default:
                                  21]
  -S, --sketch-size INTEGER       Mash parameter for Sketch size  [default:
                                  2000]
  -r, --cluster-resolution FLOAT  Resolution value used for clustering.
                                  [default: 0.75]
  -e, --eval-identity FLOAT       Blast percentage identity to use for evaluation of identification.
                                  [default: 85]
  --seed INTEGER                  Seed value for clustering
  --cpus INTEGER                  Number of cpus used for execution  [default:
                                  1]
  --force / --no-force            Overwrite the output folder if it exists
                                  [default: no-force]
  --help                          Show this message and exit.
```

Example for reference-alleles using defaults:

```bash
taranys reference-alleles \
--schema schema_dir \
--output output_reference_alleles_dir \
--eval-cluster \
--cpus number_of_cpus \
--force overwrite output dir
```

Command example changing params:

```bash
taranys reference-alleles \
--schema schema_dir \
--output output_reference_alleles_dir \
--eval-cluster \
--kmer-size k-mer size for mash \
--sketch-size Sketch size for mash \
--cluster-resolution resolution used for clustering \
--cpus number_of_cpus \
--force overwrite output dir
```

### **allele_calling:**

La llamada de alelos es la función principal de Taranis con la que se realiza la tipificación propiamente dicha. Utilizando este módulo se identifican los locus del esquema presentes en las muestras analizadas. Para ello se utilizan como query los alelos de referencia identificados anteriormente y se realiza un alineamiento con blast utilizando como base de datos los ensamblados en formato fasta. Este alineamiento nos permite obtener el alelo que está presente en la muestra, realizando una clasificación basada en las categorías descritas por el software [chewBBACA](https://chewBBACA.readthedocs.io/en/latest/user/modules/AlleleCall.html#outputs).

![allele_calling](assets/allele_calling.png)

Usage:

```bash
Usage: taranys allele-calling [OPTIONS] ASSEMBLIES...

Options:
  -s, --schema PATH               Directory where the schema with the core
                                  gene files are located.   [required]
  -r, --reference PATH            Directory where the schema reference allele
                                  files are located.   [required]
  -a, --annotation PATH           Annotation file.   [required]
  -t, --hit_lenght_perc FLOAT     Threshold value to consider in blast hit
                                  percentage regarding the reference length.
                                  Values from 0 to 1. [default:
                                  0.8]
  -p, --perc-identity INTEGER     Percentage of identity to consider in blast.
                                  [default: 85]
  -o, --output PATH               Output folder to save reference alleles
                                  [required]
  --force / --no-force            Overwrite the output folder if it exists
                                  [default: no-force]
  --snp / --no-snp                Create SNP file for alleles in assembly in
                                  relation with reference allele  [default:
                                  no-snp]
  --alignment / --no-alignment    Create alignment files  [default: no-
                                  alignment]
  -q, --proteine-threshold INTEGER
                                  Threshold of protein coverage to consider as
                                  TPR  [default: 80]
  -i, --increase-sequence INTEGER
                                  Increase the number of triplet sequences to
                                  find the stop codon  [default: 20]
  --cpus INTEGER                  Number of cpus used for execution  [default:
                                  1]
  --help                          Show this message and exit.
```
  
Command example for allele calling using defaults:

```bash
taranys allele-calling \
-s schema_dir \
-a annotation_file \
-r reference_alleles_dir \
-o output_allele_calling_dir \
--snp \
--cpus 10 \
samples_dir
```

Allele calling for blast and threshold settings:

```bash
taranys allele_calling \
-s schema_dir \
-a annotation_file \
-r reference_alleles_dir \
-o output_allele_calling_dir \
-t coverage threshold to consider in blast \
-p percentage of identity to consider in blast \
-q threshold to consider as TPR \
-i increase number of nucleotides to find stop codon \
--snp Create SNP file \
--cpus 10 \
--alignment\
samples_dir 
```

### **distance_matrix:**

The similarity between two or more genomes is estimated by comparing their respective allelic profiles and calculating the total number of differing alleles. These allelic differences between genomes are obtained by generating a distance matrix using the distance_matrix module, which takes as input the allele matrix resulting from the allele_calling process. The Hamming distance is then calculated between pairs of samples.

Usage:

```bash
Usage: taranys distance-matrix [OPTIONS]

Options:
  -a, --alleles PATH              Alleles matrix file from which to obtain
                                  distances between samples  [required]
  -o, --output PATH               Output folder to save distance matrix
                                  [required]
  --force / --no-force            Overwrite the output folder if it exists
                                  [default: no-force]
  -l, --locus-missing-threshold INTEGER
                                  Maximum percentaje of missing values a locus
                                  can have, otherwise is filtered. By default
                                  core genome is calculated, locus must be
                                  found in all samples.  [default: 0]
  -s, --sample-missing-threshold INTEGER
                                  Maximum percentaje for missing values a
                                  sample can have, otherwise it is filtered
                                  [default: 20]
  --paralog-filter / --no-paralog-filter
                                  Consider paralog tags (NIPH, NIPHEM) as
                                  missing values.  [default: paralog-filter]
  --lnf-filter / --no-lnf-filter  Consider LNF as missing values.  [default:
                                  lnf-filter]
  --plot-filter / --no-plot-filter
                                  Consider PLOT as missing values.  [default:
                                  plot-filter]
  --help                          Show this message and exit.
```

Command example:

```bash
taranys distance-matrix \
-alleles allele_calling_match.csv file \
--output distance_matrix_dir
--force overwrite output folder
```

Distance matrix with threshold settings:

```bash
taranys distance-matrix \
--alleles allele_calling_match.csv file \
--output distance_matrix_dir
--locus-missing-threshold threshold for missing locus \
--sample-missing-threshold threshold for missing samples \
--paralog-filter \
--lnf-filter \
--plot-filter \
--force overwrite output folder
```

## Output

- **analyze_schema:**

  - **FOLDERS and FILES structure:**

    - **new_schema**  Contains the new schema.
    - **prokka** Contains prokka results
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

    - **allele_annotation.tar.gz** Annotation schema file to be inputted in allele calling module.

- **reference_alleles:**

  - **FOLDERS and FILES structure:**
  
    - **Clusters** Contains the cluster allele files
      - **[cluster_alleles].txt** cluster allele file
    - **evaluate_cluster**
      - **cluster_evaluation.csv** Evaluation result with the following info:
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

- **allele_calling:**

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
    - **[locus_name]_snp_data.csv** One file per sample
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

- **distance_matrix:**

  - **FILES:**
    - **allele_matrix_fil.tsv:** Filtered allele calling matrix filtered
    - **distance_matrix_core.tsv:** distance matrix with only core genome (only locus present in all samples)
    - **distance_matrix.tsv:** distance matrix with all locus present in the schema
