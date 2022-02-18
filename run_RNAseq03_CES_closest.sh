#SBATCH --partition=genetics_1             # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=RNAseqCES          # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                  # Total # of tasks across all nodes
#SBATCH --cpus-per-task=8            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                 # Real memory (RAM) required (MB)

#SBATCH --time=96:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=RNAseqCES.%A_%a.out     # STDOUT output file
#SBATCH --error=RNAseqCES.%A_%a.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

module load bedtools2
#Sort rlog-normalized reads by DE group
bedtools intersect -wa -c -a ../RNA/RNAseq_htseq_rlog_normalized.bed -b bed/nonDE.all.bed bed/F.all.bed bed/M.all.bed | awk '$13 == "0"' | awk '$1 == "2L";$1 == "2R";$1 == "3L";$1 == "3R";$1 == "X"' | cut -f 1-4 > dmel-all-r6.27.all.expressed.Noloop.bed

bedtools intersect -wa -c -a ../RNA/RNAseq_htseq_rlog_normalized.bed -b bed/nonDE.all.bed bed/F.all.bed bed/M.all.bed | awk '$13 != "0"' | awk '$1 == "2L";$1 == "2R";$1 == "3L";$1 == "3R";$1 == "X"' | cut -f 1-4 > dmel-all-r6.27.all.expressed.nonDE.bed

bedtools intersect -wa -c -a ../RNA/RNAseq_htseq_rlog_normalized.bed -b bed/F.all.bed | awk '$13 != "0"' | awk '$1 == "2L";$1 == "2R";$1 == "3L";$1 == "3R";$1 == "X"' | cut -f 1-4 > dmel-all-r6.27.all.expressed.F.bed

bedtools intersect -wa -c -a ../RNA/RNAseq_htseq_rlog_normalized.bed -b bed/M.all.bed | awk '$13 != "0"' | awk '$1 == "2L";$1 == "2R";$1 == "3L";$1 == "3R";$1 == "X"' | cut -f 1-4 > dmel-all-r6.27.all.expressed.M.bed

#Measure distance from X genes +/- R-loops to CES
bedtools sort -i dmel-all-r6.27.all.expressed.Noloop.bed | bedtools closest -d -a - -b ~/Reference/dmel_r6_CES2.bed | awk '$1 == "X"' | sort -k 9,9n | awk {'print "Noloop\t"$4"\t"$9'} > Noloop_CES_closest.txt
bedtools sort -i dmel-all-r6.27.all.expressed.nonDE.bed | bedtools closest -d -a - -b ~/Reference/dmel_r6_CES2.bed | awk '$1 == "X"' | sort -k 9,9n | awk {'print "nonDE\t"$4"\t"$9'} > nonDE_CES_closest.txt
bedtools sort -i dmel-all-r6.27.all.expressed.F.bed | bedtools closest -d -a - -b ~/Reference/dmel_r6_CES2.bed | awk '$1 == "X"' | sort -k 9,9n | awk {'print "FE\t"$4"\t"$9'} > F_CES_closest.txt
bedtools sort -i dmel-all-r6.27.all.expressed.M.bed | bedtools closest -d -a - -b ~/Reference/dmel_r6_CES2.bed | awk '$1 == "X"' | sort -k 9,9n | awk {'print "ME\t"$4"\t"$9'} > M_CES_closest.txt

#Import into R, rbind, make boxplots
