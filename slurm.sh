#!/bin/bash

# set the number of nodes
#SBATCH --nodes=1

# cores per node
#SBATCH --ntasks-per-node=32

# set max wallclock time
#SBATCH --time=00:30:00

# set name of job
#SBATCH --job-name=ParameterEstimation

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=END,FAIL

# send mail to this address
#SBATCH --mail-user=matt.jaquiery@psy.ox.ac.uk

module purge
module load git
module load R/4.0.2-foss-2020a

export R_LIBS=$HOME/R_Libs

Rscript /home/wolf5224/arc_thesis-parameter-estimation/main.R