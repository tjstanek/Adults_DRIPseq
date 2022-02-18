#!/bin/sh

#SBATCH --partition=genetics_1       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=states          # Assign a short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=2            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                   # Real memory (RAM) required (MB)
#SBATCH --array=1-7                  # Range of jobs to run

#SBATCH --time=16:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=states.%A.%a.out  # STDOUT output file
#SBATCH --error=states.%A.%a.err   # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env


module load bedtools2

# intersect shared peaks with features
bedtools intersect -wa -wb -a  ~/Reference/Kc_chromatin_states_r6.bed -b ${SLURM_ARRAY_TASK_ID}-idr.dsA.bed | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > ${SLURM_ARRAY_TASK_ID}.idr.states.cf

# add missing features with 0 counts, sort
python append_missing_feats_states.py ${SLURM_ARRAY_TASK_ID}.idr.states.cf ${SLURM_ARRAY_TASK_ID}.idr.all.cf
sort -k1,1 -k2,2 ${SLURM_ARRAY_TASK_ID}.idr.all.cf > ${SLURM_ARRAY_TASK_ID}.idr.states.cf
rm ${SLURM_ARRAY_TASK_ID}.idr.all.cf

# create shuffled peaks 10000 times
for i in {1..10000}
do
    bedtools shuffle -chrom -g r6.22.chrom.sizes -i ${SLURM_ARRAY_TASK_ID}-idr.dsA.bed > ${SLURM_ARRAY_TASK_ID}.idr.sfl
    bedtools intersect -wa -wb -a ~/Reference/Kc_chromatin_states_r6.bed -b ${SLURM_ARRAY_TASK_ID}.idr.sfl | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k 2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > ${SLURM_ARRAY_TASK_ID}.shuffle.single.cf

    python append_missing_feats_states.py ${SLURM_ARRAY_TASK_ID}.shuffle.single.cf ${SLURM_ARRAY_TASK_ID}.shuffle.single.all.cf
    sort -k1,1 -k2,2 ${SLURM_ARRAY_TASK_ID}.shuffle.single.all.cf >> ${SLURM_ARRAY_TASK_ID}.shuffle.states.cf
    rm ${SLURM_ARRAY_TASK_ID}.shuffle.single.all.cf
done
rm ${SLURM_ARRAY_TASK_ID}.idr.sfl
rm ${SLURM_ARRAY_TASK_ID}.shuffle.single.cf

# group features, count means by feature per Autosome-vs-X
sort -k 1,2 ${SLURM_ARRAY_TASK_ID}.shuffle.states.cf | bedtools groupby -i  - -g 1-2 -c 3 -o sum | awk {'print $1"\t"$2"\t"$3/10000'} > ${SLURM_ARRAY_TASK_ID}.shuffle.idr.states.cf

# calculate pvalues for each feature
python features_pvals_states.py ${SLURM_ARRAY_TASK_ID}.shuffle.states.cf ${SLURM_ARRAY_TASK_ID}.idr.states.cf ${SLURM_ARRAY_TASK_ID}.pvals.states.cf
sort -k1,1 -k2,2 ${SLURM_ARRAY_TASK_ID}.pvals.states.cf > ${SLURM_ARRAY_TASK_ID}.pvals.idr.states.cf

rm ${SLURM_ARRAY_TASK_ID}.pvals.states.cf
