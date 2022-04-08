
# Step 1: Build the reference index
bowtie2-build --threads 33 input/reference/reference.fasta.gz tmp/reference 1> tmp/reference_indexing.log 2> tmp/reference_indexing.err

# Step 2: Map each read pair to the assembly

while read -r sample; do
# local alignment
   bowtie2 --threads 20 --local -x tmp/reference \
           -1 input/reads/${sample}.R1.fastq.gz \
           -2 input/reads/${sample}.R2.fastq.gz \
           | samtools view -b - > tmp/mappings/${sample}.bam
done <<< "$samples"

