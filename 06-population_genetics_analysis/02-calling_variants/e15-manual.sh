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

# Step 7: copy bam and bai in archive from tmp2
rsync -avx --human-readable --progress tmp2/E15-P.rg.sorted.fixmate.position.markdup.bam.bai ~/archive/2019-03-05-association_analysis_flye_assembly/2019-03-07-variant_calling/.
rsync -avx --human-readable --progress tmp/E15-P.rg.sorted.fixmate.position.markdup.bam ~/archive/2019-03-05-association_analysis_flye_assembly/2019-03-07-variant_calling/.

# Step 8: soft links to result
cd result

ln -s ~/archive/2019-03-05-association_analysis_flye_assembly/2019-03-07-variant_calling/E15-P.rg.sorted.fixmate.position.markdup.bam.bai .

ln -s ~/archive/2019-03-05-association_analysis_flye_assembly/2019-03-07-variant_calling/E15-P.rg.sorted.fixmate.position.markdup.bam .

cd ..

cd tmp

ln -s ~/archive/2019-03-05-association_analysis_flye_assembly/2019-03-07-variant_calling/E15-P.rg.sorted.fixmate.position.markdup.bam.bai .

ln -s ~/archive/2019-03-05-association_analysis_flye_assembly/2019-03-07-variant_calling/E15-P.rg.sorted.fixmate.position.markdup.bam .
