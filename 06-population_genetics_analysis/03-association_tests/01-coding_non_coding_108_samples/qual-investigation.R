library(ggplot2)
qual <- read.csv("vcf-qual3.txt", header = FALSE, sep = "")
colnames(qual) <- c("num_scaffolds", "phred_qual")
my_plot <- ggplot(data = qual, aes (x = qual$phred_qual, y = qual$num_scaffolds)) + geom_bar(stat = "identity") + theme_classic()

pdf("qual-hist7.pdf")
my_plot
dev.off()

summary(qual$phred_qual)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#     0.0       0.0       0.0     470.9       8.7 2374020.0 


## depth investigation
depth <- read.csv("vcf-depth-mean.txt", header = FALSE, sep = "")
colnames(depth) <- c("scaffold", "position", "mean_depth", "depth_variance")
depth$mean_depth <- as.numeric(as.character(depth$mean_depth))
depth$depth_variance <- as.numeric(as.character(depth$depth_variance))
short_depth <- head(depth, 2000)


pdf("depth-plot.pdf")
plot(depth$mean_depth, main="Mean Read Depth for all SNPs", xlab = "Mean depth", ylab = "Number of loci")
dev.off()

depth_plot <- ggplot(data = depth, aes (x = depth$mean_depth) + geom_bar(stat = "identity") + theme_classic()
