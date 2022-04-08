#!/bin/bash

# get clean names of all samples
samples=`cat list_of_samples_0001.txt | cut -d "_" -f 1 | sort -u`

# create file to put the results
touch coverage_per_sample.txt

# looping through samples 
while read -r sample; do
 echo ${sample} >> coverage_per_sample.txt
 grep -A 1 -E ${sample} 115_pheidole_coverage.txt | cut -f 2 | sort | head -n 2 | awk '{sum+=$1 *150 / 287000000} END {print sum}' >> coverage_per_sample.txt
done <<< "$samples"


# sum number of sequences, time it by read length, divide by genome size estimation 
cat 115_pheidole_coverage.txt | cut -f 2 | sort -n | tail -n230 | awk '{sum += $1}END{print sum*150 / 300000000}'
# 1381.47

# calculating number of sequences
cat 115_pheidole_coverage.txt | cut -f 2 | sort -n | tail -n230 | awk '{sum += $1}END{print sum}'         
#2762930432

# calculating number of nucleotides
cat 115_pheidole_coverage.txt | cut -f 2 | sort -n | tail -n230 | awk '{sum += $1}END{print sum*150}'
#414439564800
# 414439564800 / 300000000
# 1381.465


# calculating coverage
# in average, there are more sequences in M than P
grep -E -A 1 --no-group-separator "\-P" coverage_per_sample.txt | sort | grep -E -v "P" | awk '{sum+= $1} END {print sum / NR}'
# 12.2237
grep -E -A 1 --no-group-separator "\-M" coverage_per_sample.txt | sort | grep -E -v "M" | awk '{sum+= $1} END {print sum / NR}'
# 13.1155


# after removing adapters
grep -e "retained nucleotides" */*settings | cut -d ' ' -f 1,5 | ruby -pe 'gsub(/\/.*umber/, "")' > retained_read_counts.txt

