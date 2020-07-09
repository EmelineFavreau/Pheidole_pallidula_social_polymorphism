#!/usr/bin/env Rscript

# this script takes as input from Bedtools coveragebed
# in which each 5kb region has info about read depth and number of bases covered
# here a mean of read depth is calculated for each window (median normalised)
# the input is a table with 5 values, the average mean depth for each window

# set the arguments to use
args = commandArgs(trailingOnly=TRUE)

# file is named tmp/${1}/hist_5kb_${sample}_subset
my_file <- args[1]

# obtain colony name and contig name
colony_name <- gsub(x = my_file, pattern = ".*5kb_", replacement = "")
contig_name <- gsub(x = my_file, pattern = "tmp/|/hist.*", replacement = "")


print(contig_name)

# import contig read depth mean in 5kb sliding window
contig_coverage5kb <- read.table(args[1], header = FALSE, sep = "\t")

# name columns
colnames(contig_coverage5kb) <- c("contig", "start", "end", "read_depth", "nun_bases", "window_size", "proportion_window")

# create vector of window starts
window_start_vec <- unique(contig_coverage5kb$start)

# vector of means (the results will go in here)
mean_normalised_by_median <- c()

# looping through each window
for(position in 1:length(window_start_vec)){

 # subset table for this window
 contig_coverage5kb_subset <- subset(contig_coverage5kb, subset = start == window_start_vec[position])

 # calculate median for this window
 this_window_median <- median(c(rep(contig_coverage5kb_subset$read_depth, times = contig_coverage5kb_subset$nun_bases)))

 # divide each read depth by median (normalisation)
 contig_coverage5kb_subset$read_depth_median_normalised <- contig_coverage5kb_subset$read_depth / (this_window_median + 1)

 # calculate the mean for this window
 #this_window_mean <- mean(c(rep(contig_coverage5kb_subset$read_depth_median_normalised, times = contig_coverage5kb_subset$nun_bases)))
 this_window_mean <- mean(c(rep(contig_coverage5kb_subset$read_depth, times = contig_coverage5kb_subset$nun_bases)))

 # save the result (one mean per window)
 mean_normalised_by_median <- c(mean_normalised_by_median, this_window_mean)
}

# create name for output
output_name <- paste("tmp/", contig_name, "/mean_normalised_by_median", colony_name, sep = "")

# save output
write.table(mean_normalised_by_median, file = output_name, row.names = FALSE, col.names = FALSE)
