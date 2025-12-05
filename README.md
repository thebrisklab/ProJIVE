# ProJIVE: Probabilistic Joint and Individual Variance Explained

## Install instructions
```
devtools::install_github('github.com/thebrisklab/ProJIVE')
```

## Quick example using simulated data
The R script "ProJIVE_Example" generates synthetic data as described in the Simulations section of the manuscript "Probabilistic Joint and Individual Variance Explained (ProJIVE) for Data Integration." The script is written in R and exhibits implementation of two custom functions: 1) GenerateToyData is used to generate simulated data sets for which JIVE analysis is applicable; and 2) ProJIVE conducts Probabilistic JIVE. The functions that the script calls on are contained in the "R" folder. These files require the following R packages to function: ajive (https://github.com/idc9/r_jive), Matrix, ggplot2, reshape2, fields, mvtnorm, dplyr, xtable, MASS, extradistr, stringr, MCMCpack, and CJIVE.

## Replicate simulations from the manuscript, "Probabilistic Joint and Individual Variation Explained (ProJIVE) for Data Integration"
Scripts to replicate the simulation study described in "Probabilistic Joint and Individual Variance Explained (ProJIVE) for Data Integration" are located in the "Simulations" folder. Simulations were conducted on the High Power Computing cluster (HPC) at Rollins School of Public Health, Emory University. The HPC runs the CentOS Linux operating system (currently version 8.9) and uses the SLURM job scheduler. Scripts are set up to run one replication per available compute node and stores results from each replication as a CSV file. The bash scripts "create_and_submit_ProJIVE_Sims_GG.submit" and "create_and_submit_ProJIVE_Sims_MixRad.submit" build and execute bash scripts for each replicate of each setting described in the simulation study. Each script built by the "create_and_submit" scripts will submit a job that calls on either "FnRunProJIVE_Simulations_n1000_GG.R" or "FnRunProJIVE_Simulations_n1000_MixRad.R" and contains input parameters determined by the simulation settings. Lastly, the scripts "ExamineProJIVE_SimulationResults_GG.R" and "ExamineProJIVE_SimulationResults_MixRad.R" are used to generate graphical and numerical summaries of the simulation studies. All of the R scripts rely on the functions, scripts and packages listed in the section below. 

