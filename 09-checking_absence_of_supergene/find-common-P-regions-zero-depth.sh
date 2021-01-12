#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=40G
#$ -j y
#$ -l h_rt=4:00:00
#$ -o logs/$JOB_NAME/$JOB_ID/

set -euo pipefail

module load bedtools
module load bedops


# find common loci with zero depth for all P samples
bedtools intersect -sorted -a tmp/template-bam/A01-P_zero_coverage_regions.BED \
	-b tmp/*-P_zero_coverage_regions.BED \
	> tmp/P_zero_coverage_regions.BED

# resulting file to be sorted by start position 
sort-bed --max-mem 5G tmp/P_zero_coverage_regions.BED > tmp/P_zero_coverage_regions.sorted.BED


# the sorted bed file needs features to be merged
bedtools merge -i tmp/P_zero_coverage_regions.sorted.BED > tmp/P_zero_coverage_regions.merged.BED

