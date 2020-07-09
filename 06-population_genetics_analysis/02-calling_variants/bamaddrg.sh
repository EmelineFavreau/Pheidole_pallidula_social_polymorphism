#!/bin/bash
#$ -cwd                        # Run the code from the current directory
#$ -j y                        # Merge the standard output and standard error
#$ -l h_rt=4:00:00             # 4 hours per task
#$ -t 1-114                    # the number of iterations
#$ -pe smp 1                   # parallel environment, 1 core per job
#$ -l h_vmem=10G               # 10G RAM per core (10G total here)
#$ -o logs/$JOB_NAME/$JOB_ID/  # redirect the standard output from the job
#$ -tc 20 	                   # take only 20 files at once

module load samtools/1.9

set -eux

sample=$(sed -n "${SGE_TASK_ID}p" 114samples.txt)

# Step 1: Add read groups to each alignment - less than 4 minutes
~/bin/bamaddrg/bamaddrg -b input/alignments/${sample}.bam -s ${sample} > tmp/${sample}.rg.bam
