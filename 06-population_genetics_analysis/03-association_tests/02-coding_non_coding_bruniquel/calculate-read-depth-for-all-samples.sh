#!/bin/bash
module load samtools
module load R


for(( i=1; i<=114; i++))
    do
    # loop through each sample
    sample=$(sed -n "${i}p" samples-without-E15)

    mkdir -p tmp/${1}
    mkdir -p result/${1}

    #samtools view -bh input/${sample}.rg.sorted.fixmate.position.markdup.bam ${1}_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon > tmp/${1}/${sample}_subset.bam

    # obtain report for read depth by 5kb window
    # -hist report: depth of chunk | num_bases | size_of_chunk | %ofchunkatdepth
    #/data/home/btx077/software/bedtools2/bin/coverageBed -hist -sorted -a tmp/${1}_5kb.BED -b tmp/${1}/${sample}_subset.bam > tmp/$1/hist_5kb_${sample}

    # remove info about all windows
    #grep -v "all" tmp/${1}/hist_5kb_${sample} > tmp/${1}/hist_5kb_${sample}_subset

    # concatenate all actions into one liner to avoid running out of space
    samtools view -bh input/${sample}.rg.sorted.fixmate.position.markdup.bam ${1}_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon | /data/home/btx077/software/bedtools2/bin/coverageBed -hist -sorted -a tmp/${1}_5kb.BED -b stdin | grep -v "all" > tmp/${1}/hist_5kb_${sample}_subset

    # calculate read depth mean (normalised by median) for each window
    Rscript calculate-mean.R tmp/${1}/hist_5kb_${sample}_subset

    # remove temp files
    #rm -f tmp/${1}/${sample}_subset.bam
    #rm -f tmp/$1/hist_5kb_${sample}

done
