#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -j y
#$ -l h_rt=4:00:00
#$ -t 1-73
#$ -o logs/$JOB_NAME/$JOB_ID/

set -euo pipefail

sample=$(sed -n "${SGE_TASK_ID}p" input/sample-list)

module load bedtools

# obtain a bed file for each P sample with zero read depth
bedtools genomecov -bga \
	-ibam input/${sample}.rg.sorted.fixmate.position.markdup.bam \
	| awk '$4 == 0' \
	| sed -r 's/(\s+)?\S+//4' \
	| sort -k1,1 -k2,2n > tmp/${sample}_zero_coverage_regions.BED



