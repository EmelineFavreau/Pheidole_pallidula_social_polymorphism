#!/bin/bash
#$ -cwd                        # Run the code from the current directory
#$ -j y                        # Merge the standard output and standard error
#$ -l h_rt=1:00:00             # 1 task - testing if this works
#$ -t 1-115                    # the number of iterations
#$ -pe smp 8                   # parallel environment, 8 cores per job
#$ -l h_vmem=10G               # 10G RAM per core (80G total here)
#$ -o logs/$JOB_NAME/$JOB_ID/  # redirect the standard output from the job
#$ -tc 20 	                   # take only 20 files at once

sample=$(sed -n "${SGE_TASK_ID}p" samples)

# calculate read depth by position, select for two contigs
# calculate the mean by sample for those two contigs
/data/home/btx077/software/bedtools2/bin/genomeCoverageBed -d -ibam tmp/${sample}-subset.bam | grep -E "contig_1346|contig_1470" | /data/home/btx077/software/bedtools2/bin/groupBy -g 1 -c 3 -o mean > tmp/coverage-${sample}-interesting-contigs
