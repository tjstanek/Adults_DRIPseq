#!/bin/bash

#SBATCH --partition=genetics_1    # Partition (job queue)
#SBATCH --job-name=rlog        # Assign an 8-character name to your job, no spaces
#SBATCH --nodes=1                # Number of compute nodes
#SBATCH --ntasks=1               # Processes (usually = cores) on each node
#SBATCH --cpus-per-task=28       # Threads per process (or per core)
#SBATCH --export=ALL             # Export you current environment settings to the job environment
#SBATCH --time96:00:00
#SBATCH --mem=100G
#SBATCH --output=rlognorm.%N.%j.out

#rlog-normalize htseq-cts
R --no-save < RNAseq.rlog.norm.R >& RNAseq.rlog.norm.R.log

#Convert rlog-normalized output to BED format
python txt.to.bed.py

# Run bedtools intersect with various files
module load bedtools2

file=RNAseq_htseq_rlog_normalized.bed

bedtools intersect -wa -c -a $file -b ../bam/bed/nonDE.all.bed ../bam/bed/F.all.bed ../bam/bed/M.all.bed | awk '$13 == "0"' | awk '$1 == "2L";$1 == "2R";$1 == "3L";$1 == "3R";$1 == "X"' | cut -f 1-12 | awk '{print $0"\tN"}' | sed 's/2L/Autosomes/g;s/2R/Autosomes/g;s/3L/Autosomes/g;s/3R/Autosomes/g' > All_noloop.txt

bedtools intersect -wa -c -a $file -b ../bam/bed/nonDE.all.bed | awk '$13 != "0"' | awk '$1 == "2L";$1 == "2R";$1 == "3L";$1 == "3R";$1 == "X"' | cut -f 1-12 | awk '{print $0"\tY"}' | sed 's/2L/Autosomes/g;s/2R/Autosomes/g;s/3L/Autosomes/g;s/3R/Autosomes/g' > nonDE_Rloop.txt

bedtools intersect -wa -c -a $file -b ../bam/bed/F.all.bed | awk '$13 != "0"' | awk '$1 == "2L";$1 == "2R";$1 == "3L";$1 == "3R";$1 == "X"' | awk '{print $0"\tY"}' | sed 's/2L/Autosomes/g;s/2R/Autosomes/g;s/3L/Autosomes/g;s/3R/Autosomes/g' > F_Rloop.txt

bedtools intersect -wa -c -a $file -b ../bam/bed/M.all.bed | awk '$13 != "0"' | awk '$1 == "2L";$1 == "2R";$1 == "3L";$1 == "3R";$1 == "X"' | cut -f 1-12 | awk '{print $0"\tY"}' | sed 's/2L/Autosomes/g;s/2R/Autosomes/g;s/3L/Autosomes/g;s/3R/Autosomes/g' > M_Rloop.txt
