# Data Repository for "MapGL: Inferring evolutionary gain and loss of short genomic sequence features by phylogenetic maximum parsimony."

This repository can be used to reproduce source data and results contributing to all figures and tables presented in:

>*Diehl AG, and Boyle AP. MapGL: Inferring evolutionary gain and loss of short genomic sequence features by phylogenetic maximum parsimony. BMC Bioinformatics. 2020*

The primary repository for the MapGL software can be found here: https://github.com/adadiehl/mapGL.

## Instructions
Shell scripts are provided to produce each of the datasets presented. These include all steps to download and prepare ChIP-seq data, phylogenetic trees, and alignment chains prior to running MapGL. We recommend installing and running MapGL within an Anaconda environment, following instructions available within the primary git repository (https://github.com/adadiehl/mapGL).

Each shell script will retrieve all necessary data and run the appropriate MapGL command to perform a given comparison. For example, to reproduce the mammalian comparison, issue the following command after installing all dependencies:

`bash do_mammalian_analysis.sh`

Shell scripts may be run individually, in any order. Pie charts may be plotted by running commands within visualize_results.R, either within an R session, or by calling Rscript (assuming all analysis shell scripts have already been run):

`/usr/bin/env Rscript visualize_results.R`

Output will be written to the results folder.

## Dependencies
In order to reproduce the datasets and results described in the manuscript, please use the following software versions:
1. Python 3.7.6
2. MapGL 1.2.0
3. BEDTools 2.26.0
4. phast 1.3
5. R 3.6.1

## Questions and Bug Reports
Please contact Adam Diehl (adadiehl@umich.edu) with technical questions/bug reports.
