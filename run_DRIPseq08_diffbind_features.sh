#!/bin/sh

#SBATCH --partition=genetics_1       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=features          # Assign a short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=2            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                   # Real memory (RAM) required (MB)

#SBATCH --time=96:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=features.%A.%a.out     # STDOUT output file
#SBATCH --error=features.%A.%a.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env


module load bedtools2

# intersect shared peaks with features
bedtools intersect -wa -wb -a  ~/Reference/dmel_r6_combined_features_all.bed -b nonDE.all.bed | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > nonDE.diff.cf

bedtools intersect -wa -wb -a  ~/Reference/dmel_r6_combined_features_all.bed -b F.all.bed | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > F.diff.cf

bedtools intersect -wa -wb -a  ~/Reference/dmel_r6_combined_features_all.bed -b M.all.bed | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > M.diff.cf


# add missing features with 0 counts, sort
python append_missing_feats.py nonDE.diff.cf nonDE.diff.all.cf
sort -k1,1 -k2,2 nonDE.diff.all.cf > nonDE.diff.cf

python append_missing_feats.py F.diff.cf F.diff.all.cf
sort -k1,1 -k2,2 F.diff.all.cf > F.diff.cf

python append_missing_feats.py M.diff.cf M.diff.all.cf
sort -k1,1 -k2,2 M.diff.all.cf > M.diff.cf

# Create shuffled peaks 10000 times
# ...for the nonDE peakset
for i in {1..10000}
do
    bedtools shuffle -chrom -g r6.22.chrom.sizes -i nonDE.all.bed > nonDE.all.sfl
    
    bedtools intersect -wa -wb -a ~/Reference/dmel_r6_combined_features_all.bed -b nonDE.all.sfl | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > nonDE.shuffle.single.cf

    python append_missing_feats.py nonDE.shuffle.single.cf nonDE.shuffle.single.all.cf
    sort -k1,1 -k2,2 nonDE.shuffle.single.all.cf >> nonDE.shuffle.cf
    rm nonDE.shuffle.single.all.cf
done
rm nonDE.all.sfl

# ...for the Female-enriched peakset
for i in {1..10000}
do
    bedtools shuffle -chrom -g r6.22.chrom.sizes -i F.all.bed > F.all.sfl
    
    bedtools intersect -wa -wb -a ~/Reference/dmel_r6_combined_features_all.bed -b F.all.sfl | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > F.shuffle.single.cf

    python append_missing_feats.py F.shuffle.single.cf F.shuffle.single.all.cf
    sort -k1,1 -k2,2 F.shuffle.single.all.cf >> F.shuffle.cf
    rm F.shuffle.single.all.cf
done
rm F.all.sfl

# ...for the Male-enriched peakset
for i in {1..10000}
do
    bedtools shuffle -chrom -g r6.22.chrom.sizes -i M.all.bed > M.all.sfl
    
    bedtools intersect -wa -wb -a ~/Reference/dmel_r6_combined_features_all.bed -b M.all.sfl | cut -f 4-7 | sort | uniq | cut -f 1-2 | awk '$2 != "4"' | sort -k1,1 -k2,2 | uniq -c | awk {'print $2"\t"$3"\t"$1'} > M.shuffle.single.cf

    python append_missing_feats.py M.shuffle.single.cf M.shuffle.single.all.cf
    sort -k1,1 -k2,2 M.shuffle.single.all.cf >> M.shuffle.cf
    rm M.shuffle.single.all.cf
done
rm M.all.sfl

# group features, count means by feature per Autosome-vs-X
sort -k 1,2 nonDE.shuffle.cf | bedtools groupby -i - -g 1-2 -c 3 -o sum | awk {'print $1"\t"$2"\t"$3/10000'} > nonDE.shuffle.diff.cf

sort -k 1,2 F.shuffle.cf | bedtools groupby -i - -g 1-2 -c 3 -o sum | awk {'print $1"\t"$2"\t"$3/10000'} > F.shuffle.diff.cf

sort -k 1,2 M.shuffle.cf | bedtools groupby -i - -g 1-2 -c 3 -o sum | awk {'print $1"\t"$2"\t"$3/10000'} > M.shuffle.diff.cf

# calculate pvalues for each feature
python features_pvals.py nonDE.shuffle.cf nonDE.diff.cf nonDE.pvals.cf
sort -k1,1 -k2,2 nonDE.pvals.cf > nonDE.pvals.diff.cf

python features_pvals.py F.shuffle.cf F.diff.cf F.pvals.cf
sort -k1,1 -k2,2 F.pvals.cf > F.pvals.diff.cf

python features_pvals.py M.shuffle.cf M.diff.cf M.pvals.cf
sort -k1,1 -k2,2 M.pvals.cf > M.pvals.diff.cf
