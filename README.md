# ProJIVE
Probabilistic JIVE

## Simulations
Scripts to run the simulation study described in "Probabilistic Joint and Individual Variance Explained (ProJIVE) for Data Integration" are located in the "Simulations" folder. Siulations were conducted on the High Power Computing cluster (HPC) at Rollins School of Public Health, Emory University. The HPC runs the CentOS Linux operating system (currently version 8) and uses the SLURM job scheduler. Scripts are set up to run one replication per available compute node and stores results from each replication as a CSV file.

## Example using simulated data
The R script "PJIVE_ToyExample_Markdown" generates simulated data as described in the Simulations section of the manuscript "Probabilistic Joint and Individual Variance Explained (ProJIVE) for Data Integration." The script is written in RMarkdown and exhibits implementation of three custom functions: 1) GenerateToyData is used to generate simulated data sets for which JIVE analysis is applicable; 2) cc.jive conducts Canonical Joint and Individual Variation explained (currently under review for publication and available at request); and 3) ProJIVE_EM conducts Probabilistic JIVE. The functions that the script calls on are contained in "Functions_for_CJIVE" and "Functions_for_PJIVE." These files require the following R packages to function: Matrix, ggplot2, reshape2, fields, mvtnorm, dplyr, xtable, MASS, and extradistr.
