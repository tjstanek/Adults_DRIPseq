#!/bin/bash

#SBATCH --partition=genetics_1       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=deeptools          # Assign a short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=8            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                   # Real memory (RAM) required (MB)

#SBATCH --time=2:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=deeptools.%N.%j.out     # STDOUT output file
#SBATCH --error=deeptools.%N.%j.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

#ID genic R-loops by DE group 
module load bedtools2
bedtools intersect -wa -c -a bed/F.autosomes.bed -b dmel-all-r6.27.autosomes.expressed.gtf | awk '$4 != "0"' > dsA.F.autosomes.expressed.gtf
bedtools intersect -wa -c -a bed/F.X.bed -b dmel-all-r6.27.X.expressed.gtf | awk '$4 != "0"' > dsA.F.X.expressed.gtf

bedtools intersect -wa -c -a bed/M.autosomes.bed -b dmel-all-r6.27.autosomes.expressed.gtf | awk '$4 != "0"' > dsA.M.autosomes.expressed.gtf
bedtools intersect -wa -c -a bed/M.X.bed -b dmel-all-r6.27.X.expressed.gtf | awk '$4 != "0"' > dsA.M.X.expressed.gtf

bedtools intersect -wa -c -a bed/nonDE.autosomes.bed -b dmel-all-r6.27.autosomes.expressed.gtf | awk '$4 != "0"' | cut -f 1-3 > dsA.nonDE.autosomes.expressed.gtf
bedtools intersect -wa -c -a bed/nonDE.X.bed -b dmel-all-r6.27.X.expressed.gtf | awk '$4 != "0"' > dsA.nonDE.X.expressed.gtf

#Merge genic R-loops by DE group
cat dsA.F.autosomes.expressed.gtf dsA.M.autosomes.expressed.gtf dsA.nonDE.autosomes.expressed.gtf > dsA.diff.autosomes.expressed.gtf
cat dsA.F.X.expressed.gtf dsA.M.X.expressed.gtf dsA.nonDE.X.expressed.gtf > dsA.diff.X.expressed.gtf

cat dsA.F.autosomes.expressed.gtf dsA.F.X.expressed.gtf > dsA.Females.all.expressed.gtf
cat dsA.M.autosomes.expressed.gtf dsA.M.X.expressed.gtf > dsA.Males.all.expressed.gtf
cat dsA.nonDE.autosomes.expressed.gtf dsA.nonDE.X.expressed.gtf > dsA.nonDE.all.expressed.gtf

# #samtools merge replicate bams, make bigwigs
module load samtools 
for i in {3..6}
do
    samtools merge $i.merged_dsA.bam $i.r1_dsA.bam $i.r2_dsA.bam
    samtools index $i.merged_dsA.bam
done

# BAMCOMPARE
for i in {3..6}
do
    bamCompare --effectiveGenomeSize 129941135 --ignoreForNormalization X 4 --skipNonCoveredRegions -p 20 -b1 $i.merged_dsA.bam -b2 $i.Input_dsA.bam -o $i.merged.bw --binSize 20 --smoothLength 60
done

#computeMatrix-DE.expressed.gtf, plot metaprofiles
#All genic R-loops, consensus peakset, autosomes vs X
computeMatrix scale-regions --transcriptID mrna --skipZeros -p 20 -S 5.merged.bw 6.merged.bw 3.merged.bw 4.merged.bw -R dsA.diff.autosomes.expressed.gtf dsA.diff.X.expressed.gtf -b 1000 -a 1000 -o profile1.gz --regionBodyLength 3000 --samplesLabel M379 M732 F379 F732 --binSize 50 >& computeMatrix.log

plotProfile -m profile1.gz --colors black red -out vMF_DRIPseq_Adults_diff_line.png
plotHeatmap -m profile1.gz --colors black red -out vMF_DRIPseq_Adults_diff_tornado.png

#Female genic peaks
computeMatrix scale-regions --transcriptID mrna --skipZeros -p 20 -S 3.merged.bw 4.merged.bw 5.merged.bw 6.merged.bw -R dsA.Females.all.expressed.gtf -b 1000 -a 1000 -o F.profile2.gz --regionBodyLength 3000 --samplesLabel F379 F732 M379 M732 --binSize 50 >& computeMatrix.log

plotProfile -m F.profile2.gz --yMin -0.3 --yMax 0.9 -out DRIPseq_Female_merged_peaks_body_line.png
plotHeatmap -m F.profile2.gz --yMin -0.3 --yMax 0.9 --zMin -2.2 --zMax 2.3 -out DRIPseq_Female_merged_peaks_body_tornado.png --legendLocation none

#Male genic peaks
computeMatrix scale-regions --transcriptID mrna --skipZeros -p 20 -S 3.merged.bw 4.merged.bw 5.merged.bw 6.merged.bw -R dsA.Males.all.expressed.gtf -b 1000 -a 1000 -o M.profile3.gz --regionBodyLength 3000 --samplesLabel F379 F732 M379 M732 --binSize 50 >& computeMatrix.log

 plotProfile -m M.profile3.gz --yMin -0.3 --yMax 0.9 -out DRIPseq_Male_merged_peaks_body_line.png
plotHeatmap -m M.profile3.gz --yMin -0.3 --yMax 0.9 --zMin -2.2 --zMax 2.3 -out DRIPseq_Male_merged_peaks_body_tornado.png --legendLocation none

#nonDE genic peaks
computeMatrix scale-regions --transcriptID mrna --skipZeros -p 20 -S 3.merged.bw 4.merged.bw 5.merged.bw 6.merged.bw -R dsA.nonDE.all.expressed.gtf -b 1000 -a 1000 -o N.profile4.gz --regionBodyLength 3000 --samplesLabel F379 F732 M379 M732 --binSize 50 >& computeMatrix.log

plotProfile -m N.profile4.gz --yMin -0.3 --yMax 0.9 -out DRIPseq_nonDE_merged_peaks_body_line.png
plotHeatmap -m N.profile4.gz --yMin -0.3 --yMax 0.9 --zMin -2.2 --zMax 2.3 -out DRIPseq_nonDE_merged_peaks_body_tornado.png --legendLocation none
