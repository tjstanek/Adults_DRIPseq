#!/bin/bash

#SBATCH --partition=genetics_1       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=deeptools          # Assign a short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=8            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                   # Real memory (RAM) required (MB)

#SBATCH --time=4:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=deeptools.%N.%j.out     # STDOUT output file
#SBATCH --error=deeptools.%N.%j.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

multiBamSummary BED-file --genomeChunkSize 129941135 -b 1.merged_ov.bam 3.merged_dsA.bam 4.merged_dsA.bam -l Ovaries F379 F732 -o DRIPseq_multiBam_nonDE.npz --BED dsA_features/nonDE.all.bed --outRawCounts multiBam_nonDE.tab
multiBamSummary BED-file --genomeChunkSize 129941135 -b 1.merged_ov.bam 3.merged_dsA.bam 4.merged_dsA.bam -l Ovaries F379 F732 -o DRIPseq_multiBam_FE.npz --BED dsA_features/F.all.bed --outRawCounts multiBam_FE.tab
