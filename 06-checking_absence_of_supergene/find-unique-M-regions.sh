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

# output loci with M coverage but no P coverage
#bedtools intersect \
#				   -a tmp/P_zero_coverage_regions.merged.BED \
#				   -b tmp/M_non_zero_coverage_regions.merged.BED \
#				   > tmp/unique_M_regions.BED

# resulting file to be sorted by start position 
sort-bed --max-mem 5G tmp/unique_M_regions.BED > tmp/unique_M_regions.sorted.BED

# merge overlapping positions
bedtools merge -i tmp/unique_M_regions.sorted.BED > tmp/unique_M_regions.merged.BED