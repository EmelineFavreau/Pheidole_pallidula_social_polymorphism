 bwa mem -t 20 -M tmp/genome/reference.fa input/rnaseq/SRR1325085_1.fastq.gz input/rnaseq/SRR1325085_2.fastq.gz | samtools view -Sb - > tmp/map/SRR1325085.bam
 bwa mem -t 20 -M tmp/genome/reference.fa input/rnaseq/SRR1325086_1.fastq.gz input/rnaseq/SRR1325086_2.fastq.gz | samtools view -Sb - > tmp/map/SRR1325086.bam
 bwa mem -t 20 -M tmp/genome/reference.fa input/rnaseq/SRR1325087_1.fastq.gz input/rnaseq/SRR1325087_2.fastq.gz | samtools view -Sb - > tmp/map/SRR1325087.bam
 bwa mem -t 20 -M tmp/genome/reference.fa input/rnaseq/SRR1325088_1.fastq.gz input/rnaseq/SRR1325088_2.fastq.gz | samtools view -Sb - > tmp/map/SRR1325088.bam
