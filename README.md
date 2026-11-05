# taxPProValidation

## Introduction

Bugphyzz is a database of microbial annotations harmonized from different sources with a controlled vocabulary and ontology terms. Furhtermore, these annotations were propagated to uncharacterized taxa using Ancestral State Reconstruction (ASR).

This repository contains the code for reproducing the 10-fold cross-validation of the ASR methods used for propagation.

## Contents

### Main directory

Text file: `validation_summary.tsv`. This is the "final" file that 
will be reported in papers, etc. Everytime this file is updated, it most be copied to the `waldronlab/bugphyzz/inst/extdata/` repository and the `waldronlab/bugphyzz/inst/scripts/README.md` file must be updated accordingly. This must be done since the
`waldronlab/bugphyzz` package uses this file to attach validation
results to the bugphyzz annotations.

Script: `create_summary_table.R`. This script is used to generate
the `validation_summmary.tsv` file. It must be run after all of the code
in the README files of the three directories mentioned below has been run.

Directories: `discrete-binary/`, `discrete-multistate-intersection/`, and `numeric/`. These contain the outputs of the progagation processs, including the files for each test and prediction fold (1 through 10). It also contains the scripts to get these folds, calculate MCC and R2 values, logs, and a table with the results.

The README files in each of the three directories mentioned above contain instructions about how to run the code.

## Steps for reproducing

### 1. Create result files

Go inside each of the `discrete-binary/`, `discrete-multistate-intersection/`, and `numeric/` directories and follow the instructions in their respective README.md files.

Some *.csv files must have been created in each of the directories. Proceed to step 2.

### 2. Generate a summary file

In the main directory, run in a linux-like terminal:

```bash
Rscript --vanilla create_summary_table.R # R in PATH
/usr/bin/Rscript --vanilla create_summary_table.R # on supermicro
```

## Versioning

The version should be updated when generating a new `validation_summary.tsv` file.

Latest Version: v1.0.5 (in progress; see issue #2)

### Version log

Use these logs in addition to the commits. The version is important here because it is referenced in `waldonlab/bugphyzz/inst/scripts/README.md`.

### Logs

#### v1.0.5
Folds are assigned per species (strains and subspecies take their species'
fold) and, in the numeric scripts, per NCBI id rather than per tree tip, so
no taxon is predicted from its own strains or from a duplicate tip of itself
(issue #1). Re-run on the 2026-11 sources (Bergey's 2026 extraction, Madin
2020 in full, BacDive from the API) with `gc_content` added to the numeric
list (issue #2).

#### v1.0.4
Update the code based on the waldronlab/bugphyzzExports script.
Compress "Folds" into a tar.gz file for numeric, binary, and multistate.

#### v1.0.3
Remove previous attempts at performing ASR, such as using castor for
discrete attributes and using data.tree and taxonomy.

#### v1.0.2
The validation was run by setting the values accepted as TRUE to 0.8 for
binary and multistate attributes. The code and summary tables were
updated accordingly.
