
# Code Supplement: Psychopathology of worthlessness in depression

This git repository contains R code to reproduce analyses from the paper:

Harrison P, Lawrence AJ, Wang S, Liu S, Xie G, Yang X, Zahn R. The Psychopathology of Worthlessness in Depression. Front Psychiatry. 2022 May 19;13:818542. doi: 10.3389/fpsyt.2022.818542. PMID: 35664464; PMCID: PMC9160466.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9160466/

This code release was tested with R version 4.1.2, bootnet version 1.5, and NetworkComparisonTest version 2.2.1. Full information is provided in `sessionInfo.txt`.

Note: There is no random seed so results will vary slightly from run to run.

## Instructions

There are two R scripts to be run in sequence:

* 01_make_networks.R
* 02_analysis.R

Note: `01_make_networks.R` is very time consuming due to bootstrap and permutation analyses, but results are logged to files in model_objects once complete.

Required libraries are specified at the top of each script and should be installed with their dependencies prior to running the analysis.

Scripts use relative paths and so assume that R's working directory will be set to the repository folder which contains all scripts and the given directory structure.

The script option variable `niter_opt` should be set to a low number (e.g. 10) to ensure libraries are installed correctly. To replicate the paper results use `niter_opt <- 10000`.

## Scripts
The first script reads data from two csv files in the `input/` folder and writes out 
permuted/bootstrapped psychometric networks to the `model_objects/` folder.

The second script analyses the intermediate objects and writes out tables/figures to the `output/` folder.

## Input

* input/item_info.csv
* input/study1.csv
* input/study2.csv

The data files `study1.csv` and `study2.csv` contain CSV format likert data coded as integers (\{1,2,3,4,5\}), with 1 row for each subject and a column for each of the 14 SCL-90 items used in the analysis. Columns are in their original order within the SCL-90 and use variable names set out in `item_info.csv`.

## Data availability
Placeholder data is provided with this code release. Please see the paper for more information on data availability.
