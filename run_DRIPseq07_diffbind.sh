#!/bin/bash

#SBATCH --partition=genetics_1    # Partition (job queue)
#SBATCH --job-name=diffbind        # Assign an 8-character name to your job, no spaces
#SBATCH --nodes=1                # Number of compute nodes
#SBATCH --ntasks=1               # Processes (usually = cores) on each node
#SBATCH --cpus-per-task=4       # Threads per process (or per core)
#SBATCH --export=ALL             # Export you current environment settings to the job environment
#SBATCH --time=96:00:00
#SBATCH --mem=120G
#SBATCH --output=diffbind.%N.%j.out

#Run diffbind on all samples
echo "Calling Adults DE peaks with diffbind..."
R --no-save < DRIPseq06.diffbind.R >& DRIPseq06.diffbind.log
echo "Done!"
