# Aim: Make 5kb bed files for all contigs that are Hymenoptera in Ppal_E

# load contig names to create bed files from
all_hymenoptera_Ppal_E_contigs <- read.table("tmp/5kb_hymenoptera_Ppal_E_contigs", stringsAsFactors = FALSE) 

# load length of contigs
intersected_contig_length_df <- read.csv("../../../2019-03-06-minion_flye_QC/Ppal_E.contig.length", header = FALSE, sep = "", stringsAsFactors = FALSE)

# name columns
colnames(intersected_contig_length_df) <- c("contig_name", "length")

# update name of contig in length df
intersected_contig_length_df$contig_name <- gsub(x = intersected_contig_length_df$contig_name, pattern = "Ppal_E.", replacement = "")

# vector of contig to create bed file for
contig_for_bed_file_vec <- all_hymenoptera_Ppal_E_contigs$V1

# size of window
window_size <- 5000

for(position in 1:length(contig_for_bed_file_vec)){
  if(intersected_contig_length_df$length[intersected_contig_length_df$contig_name == contig_for_bed_file_vec[position]] < window_size){
    # make a bed file of just one line (start and end of contig)
    # second column is start (0-base)
    second_column <- 0

    # third column is end (1-base)
    third_column <- intersected_contig_length_df$length[intersected_contig_length_df$contig_name == contig_for_bed_file_vec[position]]

    # assemble the bed file
    bed_file <- cbind(second_column, third_column)

    # update name of contig
    contig_name_updated <- paste(contig_for_bed_file_vec[position], "_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon", sep = "")

    # add the name of contig
    bed_file <- cbind(rep(contig_name_updated, length(second_column)),
                      second_column,
                      third_column)


    # give name to bed file
    bed_file_name <- paste("tmp/", contig_for_bed_file_vec[position], "_5kb.BED", sep = "")

  } else {
    # make a bed file for each contig, 5kb window
    # second column is start (0-base)
    second_column <- seq(from = 0, to = intersected_contig_length_df$length[intersected_contig_length_df$contig_name == contig_for_bed_file_vec[position]], by = window_size)

    # third column is end (1-base)
    third_column <- c(seq(from = window_size, to = intersected_contig_length_df$length[intersected_contig_length_df$contig_name == contig_for_bed_file_vec[position]], by = window_size),
                      intersected_contig_length_df$length[intersected_contig_length_df$contig_name == contig_for_bed_file_vec[position]])

    # assemble the bed file
    bed_file <- cbind(second_column, third_column)

    # update name of contig
    contig_name_updated <- paste(contig_for_bed_file_vec[position], "_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon", sep = "")

    # add the name of contig
    bed_file <- cbind(rep(contig_name_updated, length(second_column)),
                      second_column,
                      third_column)


    # give name to bed file
    bed_file_name <- paste("tmp/", contig_for_bed_file_vec[position], "_5kb.BED", sep = "")

  }
  # save in current directory
  write.table(x = bed_file, file = bed_file_name, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

}
# the name of the contig should be : "contig_1_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon"
