#!/bin/bash
#SBATCH --job-name=r_deps
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/r_deps_%j.out
#SBATCH --error=logs/r_deps_%j.err

set -euo pipefail
mkdir -p logs

module purge
module load r/4.5.0

# IMPORTANT: use the SAME library your sbatch jobs will use
export R_LIBS_USER="$HOME/R/4.5"
mkdir -p "$R_LIBS_USER"

# Clean up any stale locks from previous failed installs
rm -rf "$R_LIBS_USER/00LOCK-"*

# Avoid running package tests during install (reduces failures/time)
R -q -e 'Sys.setenv(R_TESTS=""); 
         pkgs <- c("Rcpp","dplyr","tidyr","R6","igraph","stringr","readxl","ggplot2","viridis","scales","parallel"); 
         install.packages(pkgs, repos="https://cloud.r-project.org", Ncpus=4)'