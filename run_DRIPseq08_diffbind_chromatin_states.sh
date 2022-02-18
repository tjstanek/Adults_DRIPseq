#!/bin/sh

#SBATCH --partition=genetics_1       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=features          # Assign a short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=2            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                   # Real memory (RAM) required (MB)

#SBATCH --time=16:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=features.%A.%a.out     # STDOUT output file
#SBATCH --error=features.%A.%a.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env


module load bedtools2

# intersect shared peaks with features
bedtools intersect -wa -wb -a  ~/Reference/Kc_chromatin_states_r6.bed -b nonDE.all.bed | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > nonDE.diff.states.cf

bedtools intersect -wa -wb -a  ~/Reference/Kc_chromatin_states_r6.bed -b F.all.bed | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > F.diff.states.cf

bedtools intersect -wa -wb -a  ~/Reference/Kc_chromatin_states_r6.bed -b M.all.bed | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > M.diff.states.cf


# add missing features with 0 counts, sort
python append_missing_feats_states.py nonDE.diff.states.cf nonDE.diff.all.cf
sort -k1,1 -k2,2 nonDE.diff.all.cf > nonDE.diff.states.cf

python append_missing_feats_states.py F.diff.states.cf F.diff.all.cf
sort -k1,1 -k2,2 F.diff.all.cf > F.diff.states.cf

python append_missing_feats_states.py M.diff.states.cf M.diff.all.cf
sort -k1,1 -k2,2 M.diff.all.cf > M.diff.states.cf
rm nonDE.diff.all.cf F.diff.all.cf M.diff.all.cf

# Create shuffled peaks 10000 times
# ...for the nonDE peakset
for i in {1..10000}
do
    bedtools shuffle -chrom -g r6.22.chrom.sizes -i nonDE.all.bed > nonDE.all.sfl
    
    bedtools intersect -wa -wb -a ~/Reference/Kc_chromatin_states_r6.bed -b nonDE.all.sfl | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > nonDE.shuffle.single.cf

    python append_missing_feats_states.py nonDE.shuffle.single.cf nonDE.shuffle.single.all.cf
    sort -k1,1 -k2,2 nonDE.shuffle.single.all.cf >> nonDE.shuffle.states.cf
    rm nonDE.shuffle.single.all.cf
done
rm nonDE.all.sfl
rm nonDE.shuffle.single.cf

# ...for the Female-enriched peakset
for i in {1..10000}
do
    bedtools shuffle -chrom -g r6.22.chrom.sizes -i F.all.bed > F.all.sfl
    
    bedtools intersect -wa -wb -a ~/Reference/Kc_chromatin_states_r6.bed -b F.all.sfl | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > F.shuffle.single.cf

    python append_missing_feats_states.py F.shuffle.single.cf F.shuffle.single.all.cf
    sort -k1,1 -k2,2 F.shuffle.single.all.cf >> F.shuffle.states.cf
    rm F.shuffle.single.all.cf
done
rm F.all.sfl
rm F.shuffle.single.cf

# ...for the Male-enriched peakset
for i in {1..10000}
do
    bedtools shuffle -chrom -g r6.22.chrom.sizes -i M.all.bed > M.all.sfl
    
    bedtools intersect -wa -wb -a ~/Reference/Kc_chromatin_states_r6.bed -b M.all.sfl | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > M.shuffle.single.cf

    python append_missing_feats_states.py M.shuffle.single.cf M.shuffle.single.all.cf
    sort -k1,1 -k2,2 M.shuffle.single.all.cf >> M.shuffle.states.cf
    rm M.shuffle.single.all.cf
done
rm M.all.sfl
rm M.shuffle.single.cf

# group features, count means by feature per Autosome-vs-X
sort -k 1,2 nonDE.shuffle.states.cf | bedtools groupby -i - -g 1-2 -c 3 -o sum | awk {'print $1"\t"$2"\t"$3/10000'} > nonDE.shuffle.diff.states.cf

sort -k 1,2 F.shuffle.states.cf | bedtools groupby -i - -g 1-2 -c 3 -o sum | awk {'print $1"\t"$2"\t"$3/10000'} > F.shuffle.diff.states.cf

sort -k 1,2 M.shuffle.states.cf | bedtools groupby -i - -g 1-2 -c 3 -o sum | awk {'print $1"\t"$2"\t"$3/10000'} > M.shuffle.diff.states.cf

# calculate pvalues for each feature
python features_pvals_states.py nonDE.shuffle.states.cf nonDE.diff.states.cf nonDE.pvals.states.cf
sort -k1,1 -k2,2 nonDE.pvals.states.cf > nonDE.pvals.diff.states.cf

python features_pvals_states.py F.shuffle.states.cf F.diff.states.cf F.pvals.states.cf
sort -k1,1 -k2,2 F.pvals.states.cf > F.pvals.diff.states.cf

python features_pvals_states.py M.shuffle.states.cf M.diff.states.cf M.pvals.states.cf
sort -k1,1 -k2,2 M.pvals.states.cf > M.pvals.diff.states.cf

rm nonDE.pvals.states.cf F.pvals.states.cf M.pvals.states.cf
