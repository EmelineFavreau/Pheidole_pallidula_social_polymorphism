library(tidyverse)
library(ggplot2)

# list of r2 values for all pairwise SNPs 
# that are sig associated with social type
r2_sig_snp <- read.table("tmp/test-ld-48-only-sig-analysis-triangle-concat.ld", 
                         quote="\"", comment.char="", stringsAsFactors=FALSE)

# list of r2 values for all pairwise SNPs 
# that are not sig associated with social type (100 randomisation)
r2_random_snp <-
  read.table("tmp/test-ld-48-random-analysis-triangle-contat100.ld", 
             quote="\"", comment.char="", stringsAsFactors=FALSE)

# triangle created by PLINK includes the diagonal (1 for each row and colum pair)
sum(r2_sig_snp$V1 == 1)
summary(r2_sig_snp$V1)
sum(r2_random_snp$V1 == 1)
summary(r2_random_snp$V1)

# remove 48*n instances of 1 (diagonal)
r2_sig_snp_tidy <- r2_sig_snp$V1[r2_sig_snp$V1 != 1]
r2_random_snp_tidy <- r2_random_snp$V1[r2_random_snp$V1 != 1]

# plot histogram
hist(r2_sig_snp_tidy)
hist(r2_random_snp_tidy)

# plot a two box and scatter
data <- data.frame(
  name=c( rep("Between Loci with significant association", length(r2_sig_snp_tidy)),
          rep("Between Loci with no association", length(r2_random_snp_tidy))),
  r2=c(r2_sig_snp_tidy, r2_random_snp_tidy)
)

data %>%
  ggplot( aes(x=name, y=r2, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +

  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Pairwise Linkage Disequilibrium") +
  xlab("")

# null hypothesis: means of two populations are equal
t.test(x = r2_sig_snp_tidy,y=r2_random_snp_tidy)
# t = 86.689, df = 1127.4, p-value < 2.2e-16 alternative hypothesis: different

# measuring r2, which we expect should be higher in sig SNPs
#alternative = "greater" is the alternative that x has a larger mean than y.
t.test(x = r2_sig_snp_tidy, y = r2_random_snp_tidy, alternative = "greater")
# t = 86.689, df = 1127.4, p-value < 2.2e-16
