#!/bin/bash
#SBATCH --partition=genetics_1    # Partition (job queue)
#SBATCH --job-name=samtools        # Assign an 8-character name to your job, no spaces
#SBATCH --nodes=1                # Number of compute nodes
#SBATCH --ntasks=1               # Processes (usually = cores) on each node
#SBATCH --cpus-per-task=1       # Threads per process (or per core)
#SBATCH --export=ALL             # Export you current environment settings to the job environment
#SBATCH --time=16:00:00
#SBATCH --mem=20G
#SBATCH --output=samtools.downsample.%Nn%j.out

module load samtools

#Downsample Female reads (Autosomes and X)
samtools view -b -s 0.5 -L autosomes.bed -o 3_r1.A.bam 3.r1.mq20.sorted.bam
samtools view -b -s 0.5 -L X.bed -o 3_r1.X.bam 3.mq20.sorted.bam
samtools merge 3.r1_dsA.bam 3_r1.A.bam 3_r1.X.bam
rm 3_r1.X.bam 3_r1.A.bam

samtools view -b -s 0.5 -L autosomes.bed -o 3_r2.A.bam 3.r2.mq20.sorted.bam
samtools view -b -s 0.5 -L X.bed -o 3_r2.X.bam 3.mq20.sorted.bam
samtools merge 3.r2_dsA.bam 3_r2.A.bam 3_r2.X.bam
rm 3_r2.X.bam 3_r2.A.bam

samtools view -b -s 0.5 -L autosomes.bed -o 3_Input.A.bam 3.Input.mq20.sorted.bam
samtools view -b -s 0.5 -L X.bed -o 3_Input.X.bam 3.mq20.sorted.bam
samtools merge 3.Input_dsA.bam 3_Input.A.bam 3_Input.X.bam
rm 3_Input.X.bam 3_Input.A.bam

samtools view -b -s 0.5 -L autosomes.bed -o 4_r1.A.bam 4.r1.mq20.sorted.bam
samtools view -b -s 0.5 -L X.bed -o 4_r1.X.bam 4.mq20.sorted.bam
samtools merge 4.r1_dsA.bam 4_r1.A.bam 4_r1.X.bam
rm 4_r1.X.bam 4_r1.A.bam

samtools view -b -s 0.5 -L autosomes.bed -o 4_r2.A.bam 4.r2.mq20.sorted.bam
samtools view -b -s 0.5 -L X.bed -o 4_r2.X.bam 4.mq20.sorted.bam
samtools merge 4.r2_dsA.bam 4_r2.A.bam 4_r2.X.bam
rm 4_r2.X.bam 4_r2.A.bam

samtools view -b -s 0.5 -L autosomes.bed -o 4_Input.A.bam 4.Input.mq20.sorted.bam
samtools view -b -s 0.5 -L X.bed -o 4_Input.X.bam 4.mq20.sorted.bam
samtools merge 4.Input_dsA.bam 4_Input.A.bam 4_Input.X.bam
rm 4_Input.X.bam 4_Input.A.bam

#Downsample Male reads (autosomes only)
samtools view -b -s 0.5 -L autosomes.bed -o 5_r1.A.bam 5.r1.mq20.sorted.bam
samtools view -b -L X.bed -o 5_r1.X.bam 5.mq20.sorted.bam
samtools merge 5.r1_dsA.bam 5_r1.A.bam 5_r1.X.bam
rm 5_r1.X.bam 5_r1.A.bam

samtools view -b -s 0.5 -L autosomes.bed -o 5_r2.A.bam 5.r2.mq20.sorted.bam
samtools view -b -L X.bed -o 5_r2.X.bam 5.mq20.sorted.bam
samtools merge 5.r2_dsA.bam 5_r2.A.bam 5_r2.X.bam
rm 5_r2.X.bam 5_r2.A.bam

samtools view -b -s 0.5 -L autosomes.bed -o 5_Input.A.bam 5.Input.mq20.sorted.bam
samtools view -b -L X.bed -o 5_Input.X.bam 5.mq20.sorted.bam
samtools merge 5.Input_dsA.bam 5_Input.A.bam 5_Input.X.bam
rm 5_Input.X.bam 5_Input.A.bam

samtools view -b -s 0.5 -L autosomes.bed -o 6_r1.A.bam 6.r1.mq20.sorted.bam
samtools view -b -L X.bed -o 6_r1.X.bam 6.mq20.sorted.bam
samtools merge 6.r1_dsA.bam 6_r1.A.bam 6_r1.X.bam
rm 6_r1.X.bam 6_r1.A.bam

samtools view -b -s 0.5 -L autosomes.bed -o 6_r2.A.bam 6.r2.mq20.sorted.bam
samtools view -b -L X.bed -o 6_r2.X.bam 6.mq20.sorted.bam
samtools merge 6.r2_dsA.bam 6_r2.A.bam 6_r2.X.bam
rm 6_r2.X.bam 6_r2.A.bam

samtools view -b -s 0.5 -L autosomes.bed -o 6_Input.A.bam 6.Input.mq20.sorted.bam
samtools view -b -L X.bed -o 6_Input.X.bam 6.mq20.sorted.bam
samtools merge 6.Input_dsA.bam 6_Input.A.bam 6_Input.X.bam
rm 6_Input.X.bam 6_Input.A.bam
