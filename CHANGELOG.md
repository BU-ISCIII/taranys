# Taranys Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.0.1] - 202X-XX-XX : https://github.com/BU-ISCIII/relecov-tools/releases/tag/

Same version as before for proper publish to pip.

### Credits

Code contributions to the release:

Sara Monzon - [saramonzon](https://github.com/saramonzon)

### Modules

#### Added enhancements

#### Fixes

#### Changed

#### Removed

### Requirements

## [3.0.0] - 2025-01-02 : https://github.com/BU-ISCIII/taranys/releases/tag/3.0.0

- Code refactor to create a proper python package
- Implementation of parallel computation for execution speed.
- New implementation of reference-alleles module using leiden algorithm for allele clustering.
- New param defaults based on empirical testing.
- Distance matrix is rewritten with more params available.

## [2.0.1beta] - 2021-07-14 : https://github.com/BU-ISCIII/taranys/releases/tag/2.0.1

### New features

- New default missing values threshold imposed for samples in distance matrix creation
- New default perc_identity_ref value for loci search in allele calling analysis
- Github actions for testing and docker generation.

### Bug fixes

- multiple statistic modes bug fixed
- BLAST database creation bug due to punctuation signs in file name fixed
- ST profile identification bug fixed
- Python modules installed from conda-forge in environment.yml

### Documentation

- Channel order in conda installation command

## [2.0.0beta] - 2021-06-24 : https://github.com/BU-ISCIII/taranys/releases/tag/2.0.0edge

Pre-release

## [0.3.3beta] - 2018-11-05 : https://github.com/BU-ISCIII/taranys/releases/tag/0.3.3

**BUG FIX:**

- Fix num_threads in blast commands

## [0.3.2beta] - 2018-11-02 : https://github.com/BU-ISCIII/taranys/releases/tag/0.3.2

**BUG FIX:**

- Added ERROR posibility to allele classification, when the allele is not PLOT but is too near to contig end, and protein codification finishes without finding a stop codon.

## [0.3.1beta] - 2018-10-27 : https://github.com/BU-ISCIII/taranys/releases/tag/0.3.1

- Added cpus as parameter.

## [0.3.0beta] - 2018-10-25 : https://github.com/BU-ISCIII/taranys/releases/tag/0.3.0

**Bug fixes:**

- Allow dependency check for versions greater than needed.

**Features:**

- Added graphics for schema evaluation(beta)
- Added SNP calling (beta)

## [0.0.1beta] - 2018-10-22 : https://github.com/BU-ISCIII/taranys/releases/tag/0.1

- cg/wgMLST analysis using assembled contigs as input using a defined schema. Comparison matrix is generated for phyloviz visualization. Blast hit are classified as Exact match , new allele, etc.
- cg/wgMLST statistics.
- beta: SNP analysis.
- Insertions, deletions and paralogues detection.