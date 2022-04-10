#!/bin/bash
#$ -cwd                        # Run the code from the current directory
#$ -j y                        # Merge the standard output and standard error
#$ -l h_rt=1:00:00             # 1 task - testing if this works
#$ -t 1-114                    # the number of iterations
#$ -pe smp 24                  # parallel environment, 24 cores per job
#$ -l h_vmem=10G               # 10G RAM per core (240G total here)
#$ -o logs/$JOB_NAME/$JOB_ID/  # redirect the standard output from the job
#$ -tc 20 	                   # take only 20 files at once


module load samtools/1.9

set -eux

sample=$(sed -n "${SGE_TASK_ID}p" 114samples.txt)

# Step 1: sort alignments by names - 7 minutes multithreaded
samtools sort -n -@ ${NSLOTS} tmp/${sample}.rg.bam \
              > tmp/${sample}.rg.sorted.bam

# Step 2: add ms and MC tags for markdup - 50 minutes not multithreaded
samtools fixmate -@ ${NSLOTS} -m tmp/${sample}.rg.sorted.bam \
              tmp/${sample}.rg.sorted.fixmate.bam

# Step 3: Markdup needs position order - 5 minutes multithreaded
samtools sort -@ ${NSLOTS} -o tmp/${sample}.rg.sorted.fixmate.position.bam \
              tmp/${sample}.rg.sorted.fixmate.bam

# Step 4: mark duplicates - not multithreaded
samtools markdup -@ ${NSLOTS} tmp/${sample}.rg.sorted.fixmate.position.bam \
                 tmp/${sample}.rg.sorted.fixmate.position.markdup.bam

# Step 5: obtain bam index
samtools index -@ ${NSLOTS} tmp/${sample}.rg.sorted.fixmate.position.markdup.bam \
                tmp/${sample}.rg.sorted.fixmate.position.markdup.bam.bai
