#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -j y
#$ -l h_rt=4:00:00
#$ -o logs/$JOB_NAME/$JOB_ID/

set -euo pipefail

module load bedtools

# sort bed by chromosome and by feature size (in descending order)
sortBed -i tmp/unique_M_regions.merged.BED -sizeD > tmp/longest_unique_M_regions.merged.BED

# make a 4th column with feature length
awk 'BEGIN { OFS = "\t" } { $4 = $3 - $2 } 1' \
	tmp/longest_unique_M_regions.merged.BED > tmp/contig-start-end-length.txt

# sort by feature length, key starting and ending with column 4, descending order
sort --key 4,4 --reverse --numeric-sort tmp/contig-start-end-length.txt > tmp/contig-start-end-length.sorted.txt
