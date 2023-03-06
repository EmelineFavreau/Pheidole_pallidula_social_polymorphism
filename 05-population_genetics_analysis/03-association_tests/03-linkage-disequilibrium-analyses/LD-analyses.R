#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# counter
i <- args[1]

# list of all SNPs for their Chromosome and Position
allsnpChromPos <- read.delim("tmp/allsnpChromPos", 
                             header=FALSE, 
                             stringsAsFactors=FALSE)

# name columns
colnames(allsnpChromPos) <- c("CHROM", "POS")

# list of 48 SNPs for their Chromosome and Position
chrom_pos_fortyeightsnps <- read.delim("input/chrom_pos_fortyeightsnps.txt",
                                       header=FALSE,
                                       stringsAsFactors=FALSE)

# name columns
colnames(chrom_pos_fortyeightsnps) <- c("CHROM", "POS")

# remove from all SNPs the 48 SNPs known to be associated with social supergene
tmp_df <- rbind(allsnpChromPos, chrom_pos_fortyeightsnps)

allsnpChromPosno48 <- tmp_df %>%
  group_by(across(everything())) %>%
  filter(n() == 1)

# select 48 random rows
test <- as.data.frame(allsnpChromPosno48) %>% slice_sample(n = 48)

# make a file name
filename <- paste("tmp/random48snps", i, sep = "")
# save this file
write.table(test, file = filename, quote = FALSE,
            row.names = FALSE, col.names = FALSE)
