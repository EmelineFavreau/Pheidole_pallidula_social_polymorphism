#!/bin/bash
#$ -cwd                        # Run the code from the current directory
#$ -j y                        # Merge the standard output and standard error
#$ -l h_rt=1:00:00             # 1 task - testing if this works
#$ -t 1-115                    # the number of iterations
#$ -pe smp 8                   # parallel environment, 8 cores per job
#$ -l h_vmem=5G               # 10G RAM per core (80G total here)
#$ -o logs/$JOB_NAME/$JOB_ID/  # redirect the standard output from the job
#$ -tc 20 	                   # take only 20 files at once

sample=$(sed -n "${SGE_TASK_ID}p" samples)

# prepare the output with the name of colony
echo ${sample} > tmp/${sample}_name

# calculate read depth by position, select for two contigs
/data/home/btx077/software/bedtools2/bin/genomeCoverageBed -d -ibam tmp/${sample}-subset.bam | grep -E "contig_1346|contig_1470" > tmp/read_depth_by_position_${sample}

# calculate the median by sample for those two contigs
/data/home/btx077/software/bedtools2/bin/groupBy -i tmp/read_depth_by_position_${sample} -g 1 -c 3 -o median > tmp/median_${sample}

# find the median for extreme contig
extrememedian="$(grep -E "contig_1346" tmp/median_${sample} | cut -f 2 )"

# divide each read count by median, followed by calculating the mean (output: 1 value)
awk -v m="extrememedian" '{ print $3 / (m + 1) }' tmp/read_depth_by_position_${sample}| awk '{sum+=$1} END { print sum / NR }' | paste tmp/${sample}_name - > tmp/contig_1346_mean_of_read_depth_normalised_by_median_${sample}


# find the median for normal contig
normalmedian="$(grep -E "contig_1470" tmp/median_${sample} | cut -f 2 )"

# divide each read count by median, followed by calculating the mean (output: 1 value)
awk -v m="normalmedian" '{ print $3 / (m + 1) }' tmp/read_depth_by_position_${sample}| awk '{sum+=$1} END { print sum / NR }' | paste tmp/${sample}_name - > tmp/contig_1470_mean_of_read_depth_normalised_by_median_${sample}

# for each sample, there should be two files, one for each contig
