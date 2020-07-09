#!/bin/env bash
# make a variable with all sample names
samples=`ls input/*.gz | cut -d '/' -f 2| cut -d '.' -f 1 | sort | uniq`
# looping through samples
while read -r sample; do
# extract top matches = EF518381.1 = Pheidole pallidula?
   ./run_magic_blast_em.sh input/reference input/${sample}.R1.fastq.gz input/${sample}.R2.fastq.gz
done <<< "$samples"
