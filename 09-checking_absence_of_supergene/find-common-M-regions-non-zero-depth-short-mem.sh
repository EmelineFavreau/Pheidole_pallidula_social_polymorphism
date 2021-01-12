#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -j y
#$ -l h_rt=4:00:00
#$ -o logs/$JOB_NAME/$JOB_ID/

set -euo pipefail

module load bedtools
module load bedops


# find common loci with non-zero depth for all M samples
bedtools intersect -sorted -a tmp/template-bam/A07-M_non_zero_coverage_regions.BED \
	-b tmp/*_non_zero_coverage_regions.BED \
	> tmp/M_non_zero_coverage_regions.BED

# resulting file to be sorted by start position 
sort-bed --max-mem 5G tmp/M_non_zero_coverage_regions.BED > tmp/M_non_zero_coverage_regions.sorted.BED

# merge overlapping positions
bedtools merge -i tmp/M_non_zero_coverage_regions.sorted.BED > tmp/M_non_zero_coverage_regions.merged.BED

