module load bwa/0.7.17
module load samtools/1.9
 
#######
# map pairs of reads with bwa  (STAR had no pairs on different scaffolds)
mkdir -p tmp/genome
cp input/reference.fa tmp/genome
bwa index tmp/genome/reference.fa

samples=`ls input/rnaseq | cut -f 1 -d '_' | sort | uniq`
mkdir tmp/map
echo $samples | parallel echo "bwa mem -t 20 -M tmp/genome/reference.fa input/rnaseq/{}_1.fastq.gz input/rnaseq/{}_2.fastq.gz | samtools view -Sb - > tmp/map/{}.bam" >> commands.sh
sh commands.sh > tmp/map/log 2> tmp/map/err

# they are sorted by read name; keeping this when merging (Agouti takes 1 file). Removing multipmapping reads: qual minimum 50 and those with SA:Z
samtools merge -n - tmp/map/SRR1325085.bam tmp/map/SRR1325086.bam tmp/map/SRR1325087.bam tmp/map/SRR1325088.bam | samtools view -b -q 50 - > tmp/map_merged.bam

# removing all multiple-mapping reads (needed by Agouti), and low quality mapping: 
#  UNFORTUNATELY, SOME MULTIPLE MAPPED REMAIN... UNSURE HOW TO REMOVE THEM ! (SA:Z indicates something but I'm not sure its this!)

mv tmp/map_merged.bam results/FourRnaseq_bwa_pheidoleE.bam
#####
# Then run agouti. 
#module load intelpython/2.7.12
#mkdir tmp/agouti
#python ../../soft/AGOUTI/agouti.py scaffold -assembly input/reference.fa -bam tmp/map_merged_q100.bam -outdir tmp/agouti
