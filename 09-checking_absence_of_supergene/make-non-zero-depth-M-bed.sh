#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -j y
#$ -l h_rt=4:00:00
#$ -t 1-38

sample=$(sed -n "${SGE_TASK_ID}p" input/M-sample-list)

module load bedtools

# obtain a bed file for each M sample with non-zero read depth
bedtools genomecov -bga \
	-ibam input/${sample}.rg.sorted.fixmate.position.markdup.bam \
	| awk '$4 > 0' \
	| sed -r 's/(\s+)?\S+//4' \
	| sort -k1,1 -k2,2n > tmp/${sample}_non_zero_coverage_regions.BED



