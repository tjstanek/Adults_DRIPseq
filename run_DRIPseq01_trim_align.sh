#!/bin/bash

#SBATCH --partition=genetics_1       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=dripseq          # Assign a short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=8            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                   # Real memory (RAM) required (MB)
#SBATCH --array=3-6                  # Array range
#SBATCH --time=96:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=dripseq.%A_%a.out   # STDOUT output file
#SBATCH --error=dripseq.%A_%a.err    # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

#Trim reads with trimmomatic
echo "Trimming reads..."
module load java
java -jar /home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog ${SLURM_ARRAY_TASK_ID}.r1_trim.log ${SLURM_ARRAY_TASK_ID}.r1_1.fq.gz ${SLURM_ARRAY_TASK_ID}.r1_2.fq.gz ${SLURM_ARRAY_TASK_ID}.r1_1.paired.fastq ${SLURM_ARRAY_TASK_ID}.r1_1.unpaired.fastq ${SLURM_ARRAY_TASK_ID}.r1_2.paired.fastq ${SLURM_ARRAY_TASK_ID}.r1_2.unpaired.fastq  ILLUMINACLIP:/home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:15 TRAILING:15 SLIDINGWINDOW:3:10 MINLEN:36

java -jar /home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog ${SLURM_ARRAY_TASK_ID}.r2_trim.log ${SLURM_ARRAY_TASK_ID}.r2_1.fq.gz ${SLURM_ARRAY_TASK_ID}.r2_2.fq.gz ${SLURM_ARRAY_TASK_ID}.r2_1.paired.fastq ${SLURM_ARRAY_TASK_ID}.r2_1.unpaired.fastq ${SLURM_ARRAY_TASK_ID}.r2_2.paired.fastq ${SLURM_ARRAY_TASK_ID}.r2_2.unpaired.fastq  ILLUMINACLIP:/home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:15 TRAILING:15 SLIDINGWINDOW:3:10 MINLEN:36

java -jar /home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog ${SLURM_ARRAY_TASK_ID}.Input_trim.log ${SLURM_ARRAY_TASK_ID}.Input_1.fq.gz ${SLURM_ARRAY_TASK_ID}.Input_2.fq.gz ${SLURM_ARRAY_TASK_ID}.Input_1.paired.fastq ${SLURM_ARRAY_TASK_ID}.Input_1.unpaired.fastq ${SLURM_ARRAY_TASK_ID}.Input_2.paired.fastq ${SLURM_ARRAY_TASK_ID}.Input_2.unpaired.fastq  ILLUMINACLIP:/home/ts928/Genomics_Workshop/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:15 TRAILING:15 SLIDINGWINDOW:3:10 MINLEN:36

#Align trimmed reads with bowtie
echo "Aligning reads..."
module load bowtie2
bowtie2 --no-mixed --no-discordant --dovetail --phred33 -X 1000 -q -x /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta -1 ${SLURM_ARRAY_TASK_ID}.r1_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.r1_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.r1.sam
bowtie2 --no-mixed --no-discordant --dovetail --phred33 -X 1000 -q -x /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta -1 ${SLURM_ARRAY_TASK_ID}.r2_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.r2_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.r2.sam
bowtie2 --no-mixed --no-discordant --dovetail --phred33 -X 1000 -q -x /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta -1 ${SLURM_ARRAY_TASK_ID}.Input_1.paired.fastq -2 ${SLURM_ARRAY_TASK_ID}.Input_2.paired.fastq -S ${SLURM_ARRAY_TASK_ID}.Input.sam

#Convert sam to bam
echo "Converting sam to bam..."
module load samtools
samtools view -S -b ${SLURM_ARRAY_TASK_ID}.r1.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.r1.mq20.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.r1.mq20.sorted.bam

samtools view -S -b ${SLURM_ARRAY_TASK_ID}.r2.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.r2.mq20.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.r2.mq20.sorted.bam

samtools view -S -b ${SLURM_ARRAY_TASK_ID}.Input.sam | samtools view -q 20 -b | samtools sort > ${SLURM_ARRAY_TASK_ID}.Input.mq20.sorted.bam
samtools index ${SLURM_ARRAY_TASK_ID}.Input.mq20.sorted.bam
