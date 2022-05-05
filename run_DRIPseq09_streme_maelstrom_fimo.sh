#!/bin/bash

#SBATCH --partition=genetics_1             # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=ame          # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                  # Total # of tasks across all nodes
#SBATCH --cpus-per-task=2            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120G                 # Real memory (RAM) required (MB)

#SBATCH --time=12:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=streme.%N.%j.out     # STDOUT output file
#SBATCH --error=streme.%N.%j.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

#ID motifs in each DE group using STREME
streme --p nonDE.all.fasta --oc nonDE_1e-10 --maxw 30 --pvt 1e-10
streme --p F.all.fasta --oc F_1e-2 --maxw 30 --pvt 1e-2
streme --p M.all.fasta --oc M_1e-2 --maxw 30 --pvt 1e-2

#Compile top 5 STREME motifs from each group into one PFM file
head -n 188 nonDE_1e-10/streme.txt | tail -n 159 > gimme.streme.pfm
head -n 138 F_1e-2/streme.txt | tail -n 109 >> gimme.streme.pfm
head -n 152 M_1e-2/streme.txt | tail -n 123 >> gimme.streme.pfm

#BED to gimme maelstrom input format
echo "loc" \t "cluster" > Diff.all.txt
awk {'print $1 "\t" $3'} Diff.all.txt > Diff.sex.all.txt
cp Diff.sex.all.txt Diff.strain.all.txt

awk {'print "chr" $1 ":" $2 "-" $3 "\tnonDE"'} nonDE.all.bed >> Diff.sex.all.txt
awk {'print "chr" $1 ":" $2 "-" $3 "\tFE"'} F.all.bed >> Diff.sex.all.txt
awk {'print "chr" $1 ":" $2 "-" $3 "\tME"'} M.all.bed >> Diff.sex.all.txt


#Run maelstrom on Diffbind by sex
gimme maelstrom -p gimme.streme.pfm -N 1 --no-filter Diff.sex.all.txt dm6 maelstrom.Diffsex.out

#Run FIMO to find STREME-identified motifs throughout genome
#nonDE-enriched motifs
echo "Mapping nonDE-enriched motifs..."
fimo --text --motif 1-CAAAAAAAAAAAAAAAAAAAAAAAAAAAAA nonDE_1e-10/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' > fimo/N_fimo.bed

fimo --text --motif 2-RGAGRRRGAGRRRGAGRRRGAGRRRGAGR nonDE_1e-10/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' >> fimo/N_fimo.bed

fimo --text --motif 3-CACACACACACACACACACACACACACA nonDE_1e-10/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' >> fimo/N_fimo.bed

fimo --text --motif 4-CARCARCARCARCARCARCARCARCAR nonDE_1e-10/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' >> fimo/N_fimo.bed

fimo --text --motif 5-CCMMCCAMCCMMCCMMCCMMCCMMCCMMCC nonDE_1e-10/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' >> fimo/N_fimo.bed

#Female-enriched motifs
echo "Mapping Female-enriched motifs..."
fimo --text --motif 1-CTCCTCCTCCTCCTCCTCCTCCT F_1e-2/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' > fimo/F_fimo.bed

fimo --text --motif 2-GGCGAAGGAG F_1e-2/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' >> fimo/F_fimo.bed

fimo --text --motif 3-CCACATCCACATCCACATCCAH F_1e-2/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' >> fimo/F_fimo.bed

fimo --text --motif 4-AGAGAGAGAGAGAGAGRGAGASRG F_1e-2/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' >> fimo/F_fimo.bed

fimo --text --motif 5-CAACAACAACAACAA F_1e-2/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' >> fimo/F_fimo.bed

#Male-enriched motifs
echo "Mapping Male-enriched motifs..."
fimo --text --motif 1-ACAAAAAAAAAAAAAAAAAAAAAAAAA M_1e-2/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' > fimo/M_fimo.bed

fimo --text --motif 2-ACACACACACACACACACACAC M_1e-2/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' >> fimo/M_fimo.bed

fimo --text --motif 3-AGRGAGAGRRRGAGRGRKRGRGVGAGR M_1e-2/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' >> fimo/M_fimo.bed

fimo --text --motif 4-AABGAAATGAAAWBR M_1e-2/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' >> fimo/M_fimo.bed

fimo --text --motif 5-AAAAAAAAAAAAAAAAA M_1e-2/streme.txt /projects/genetics/ellison_lab/genomes/dna/dmel-all-chromosome-r6.22.fasta | grep -v motif | awk '( NR>1) {print $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$8}' | bedtools sort -i - | bedtools merge -i - -c 4,5,6  -o distinct,mean,mean | awk '$1=="2R";$1=="2L";$1=="3R";$1=="3L";$1=="4";$1=="X"' >> fimo/M_fimo.bed
