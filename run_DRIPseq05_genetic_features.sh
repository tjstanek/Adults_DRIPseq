#!/bin/sh

#SBATCH --partition=genetics_1       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=features          # Assign a short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=2            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                   # Real memory (RAM) required (MB)
#SBATCH --array=1-7                  # Range of jobs to run

#SBATCH --time=96:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=features.%A.%a.out  # STDOUT output file
#SBATCH --error=features.%A.%a.err   # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env


module load bedtools2

# intersect shared peaks with features
bedtools intersect -wa -wb -a  ~/Reference/dmel_r6_combined_features_all.bed -b ${SLURM_ARRAY_TASK_ID}-idr.dsA.bed | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > ${SLURM_ARRAY_TASK_ID}.idr.allchrom.cf

# add missing features with 0 counts, sort
python append_missing_feats_allchrom.py ${SLURM_ARRAY_TASK_ID}.idr.allchrom.cf ${SLURM_ARRAY_TASK_ID}.idr.all.cf
sort -k1,1 -k2,2 ${SLURM_ARRAY_TASK_ID}.idr.all.cf > ${SLURM_ARRAY_TASK_ID}.idr.allchrom.cf
rm ${SLURM_ARRAY_TASK_ID}.idr.all.cf

# create shuffled peaks 10000 times
for i in {1..10000}
do
    bedtools shuffle -chrom -g r6.22.chrom.sizes -i ${SLURM_ARRAY_TASK_ID}-idr.dsA.bed > ${SLURM_ARRAY_TASK_ID}.idr.sfl
    bedtools intersect -wa -wb -a ~/Reference/dmel_r6_combined_features_all.bed -b ${SLURM_ARRAY_TASK_ID}.idr.sfl | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k 2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > ${SLURM_ARRAY_TASK_ID}.shuffle.single.cf

    python append_missing_feats_allchrom.py ${SLURM_ARRAY_TASK_ID}.shuffle.single.cf ${SLURM_ARRAY_TASK_ID}.shuffle.single.all.cf
    sort -k1,1 -k2,2 ${SLURM_ARRAY_TASK_ID}.shuffle.single.all.cf >> ${SLURM_ARRAY_TASK_ID}.shuffle.allchrom.cf
    rm ${SLURM_ARRAY_TASK_ID}.shuffle.single.all.cf
done
rm ${SLURM_ARRAY_TASK_ID}.idr.sfl
rm ${SLURM_ARRAY_TASK_ID}.shuffle.single.cf

# group features, count means by feature per Autosome-vs-X
sort -k 1,2 ${SLURM_ARRAY_TASK_ID}.shuffle.allchrom.cf | bedtools groupby -i  - -g 1-2 -c 3 -o sum | awk {'print $1"\t"$2"\t"$3/10000'} > ${SLURM_ARRAY_TASK_ID}.shuffle.idr.allchrom.cf

# calculate pvalues for each feature
python features_pvals_allchrom.py ${SLURM_ARRAY_TASK_ID}.shuffle.allchrom.cf ${SLURM_ARRAY_TASK_ID}.idr.allchrom.cf ${SLURM_ARRAY_TASK_ID}.pvals.allchrom.cf
sort -k1,1 -k2,2 ${SLURM_ARRAY_TASK_ID}.pvals.allchrom.cf > ${SLURM_ARRAY_TASK_ID}.pvals.idr.allchrom.cf

rm ${SLURM_ARRAY_TASK_ID}.pvals.allchrom.cf
