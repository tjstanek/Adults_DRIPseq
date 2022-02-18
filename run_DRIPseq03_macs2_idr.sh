#SBATCH --partition=genetics_1       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=macs2               # Assign a short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=4            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=100G                   # Real memory (RAM) required (MB)
#SBATCH --array=3-6                 # Array range

#SBATCH --time=24:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=macs2.%A_%a.out     # STDOUT output file
#SBATCH --error=macs2.%A_%a.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

macs2 callpeak --outdir idr -t ${SLURM_ARRAY_TASK_ID}.r1_dsA.bam -c ${SLURM_ARRAY_TASK_ID}.Input_dsA.bam -f BAMPE -g dm -n ${SLURM_ARRAY_TASK_ID}.r1_dsA -B -p 1e-3 2> idr/macs/${SLURM_ARRAY_TASK_ID}.r1_dsA_macs2.log
sort -k8,8nr idr/${SLURM_ARRAY_TASK_ID}.r1_dsA_peaks.narrowPeak > idr/macs/${SLURM_ARRAY_TASK_ID}.r1_dsA_peaks.narrowPeak

macs2 callpeak --outdir idr -t ${SLURM_ARRAY_TASK_ID}.r2_dsA.bam -c ${SLURM_ARRAY_TASK_ID}.Input_dsA.bam -f BAMPE -g dm -n ${SLURM_ARRAY_TASK_ID}.r2_dsA -B -p 1e-3 2> idr/macs/${SLURM_ARRAY_TASK_ID}.r2_dsA_macs2.log
sort -k8,8nr idr/${SLURM_ARRAY_TASK_ID}.r2_dsA_peaks.narrowPeak > idr/macs/${SLURM_ARRAY_TASK_ID}.r2_dsA_peaks.narrowPeak

idr --samples idr/${SLURM_ARRAY_TASK_ID}.r1_dsA_peaks.narrowPeak idr/${SLURM_ARRAY_TASK_ID}.r2_dsA_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${SLURM_ARRAY_TASK_ID}.ds-idr --plot --log-output-file ${SLURM_ARRAY_TASK_ID}.ds.idr.log

#Convert narrowPeak files into bed for diffbind
echo "Creating BED files..."
cut -f 1-6 ${SLURM_ARRAY_TASK_ID}.ds-idr > macs2/${SLURM_ARRAY_TASK_ID}-idr.dsA.bed
echo "Done creating BED files!"
