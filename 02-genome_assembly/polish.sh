#!/bin/sh
# Copyright 2017 Anurag Priyam, Queen Mary University of London.

# Load tools.
module load bwa/0.7.17
module load samtools/1.4.1

# Load common code that helps us run jobs.
source ./scripts/template.sh

# Path to the assembly file.
ASM=${params[0]}
R1=${params[1]}
R2=${params[2]}
FIX=${params[3]}
ALN=pilon/bwa_mem.bam

# Map Illumina reads to arrow polished assembly.
bwa index ${ASM}
bwa mem -t ${CPU} ${ASM} ${R1} ${R2} | samtools view -b \
  | samtools sort -@ ${CPU} -T ${TMPDIR} > ${ALN}

# Index the BAM file.
samtools index -@ ${CPU} ${ALN}

# Run pilon.
java -Xmx${RAM} -jar ${HOME}/tools/pilon-1.22.jar --threads ${CPU} \
  --genome ${ASM} --frags ${ALN} --fix ${FIX} --diploid  \
  --changes --outdir pilon > pilon/log
