#!/bin/sh
#SBATCH --job-name=ABCD_Bootstrapped_EMM
#SBATCH --mail-type=START,FAIL
#SBATCH --mail-user=redacted

#SBATCH --account=redacted
#SBATCH --partition=redacted
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00 #Max runtime in DD-HH:MM:SS format.
#SBATCH --error=/logfiles/ABCD_Bootstrapped_EMM.e%A.%a
#SBATCH --output=/logfiles/ABCD_Bootstrapped_EMM.o%A.%a # where STDOUT goes

Rscript Bootstrapped_EMM.R
