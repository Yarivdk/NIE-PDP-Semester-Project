#!/bin/bash

# Setting the control parameters of the Slurm scheduler
#SBATCH --job-name=checkpoint3
#SBATCH --nodes=1

#SBATCH --output="%x-%J.out"
#SBATCH --error="%x-%J.err"

# Retrieve assigned CPUs
NUM_CPUS=${SLURM_CPUS_PER_TASK}

# Get additional command-line arguments
ARG1=$1
ARG2=$2

# Activation of HPE CPE
source /etc/profile.d/zz-cray-pe.sh

# Setting environment variables for the scheduled task
#module load cray-mvapich2_pmix_nogpu/2.3.7
module load cray-mvapich2_pmix_nogpu

# After the srun command, write the path to your program and its arguments for running on the scheduled computing nodes:
srun ./output ${ARG1} ${ARG2} ${NUM_CPUS}

exit 0