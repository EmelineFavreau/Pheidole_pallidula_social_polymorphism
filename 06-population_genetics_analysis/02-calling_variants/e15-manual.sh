module load samtools/1.9

# Step 1: Add read groups to each alignment
~/bin/bamaddrg/bamaddrg -b input/alignments/E15-P.bam \
                        -s E15-P > tmp2/E15-P.rg.bam

# Step 2: sort alignments by names
samtools sort -n -@ 24 tmp2/E15-P.rg.bam \
              > tmp2/E15-P.rg.sorted.bam

# Step 3: add ms and MC tags for markdup to use later
samtools fixmate -m tmp2/E15-P.rg.sorted.bam \
              tmp2/E15-P.rg.sorted.fixmate.bam

# Step 4: Markdup needs position order
samtools sort -@ 24 -o tmp2/E15-P.rg.sorted.fixmate.position.bam \
              tmp2/E15-P.rg.sorted.fixmate.bam

# Step 5: mark duplicates
samtools markdup tmp2/E15-P.rg.sorted.fixmate.position.bam \
                 tmp2/E15-P.rg.sorted.fixmate.position.markdup.bam

# Step 6: obtain bam index
samtools index -@ 24 tmp2/E15-P.rg.sorted.fixmate.position.markdup.bam \
                tmp2/E15-P.rg.sorted.fixmate.position.markdup.bam.bai

