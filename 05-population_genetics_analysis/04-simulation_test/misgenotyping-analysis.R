# Emeline Favreau
# QMUL
# Pheidole project


### Misgenotyping

# Misgenotyping could have happen because of null alleles during microsatellite 
# PCR steps. We tested the power of our analysis (Fisher's exact test) for 
# isgenotyping the social type (multiple-queen or single-queen).
# 
# We simulated 1000 times our association analysis, for each time 5% of 
# the samples were assigned the alternative social type (e.g. single-queen 
# sample becomes multiple-queen sample). 


# load libraries
# install.packages("ggplot2")
# install.packages("tidyverse")
library("ggplot2")
library("tidyverse")



# import misgenotyping data

# path to the data
data_path <- "tmp"   

# list the files
files <- dir(data_path, 
             pattern = "2021-01-13-misgenotyping-5-Fisher-Bonferroni.*.assoc.fisher") 

# load data in a dataframe
data <- data_frame(filename = files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read.csv(file.path(data_path, .), sep = "", header = TRUE)) # a new data column
  ) 




# also, load the real P values to compare distributions
# Evaluating raw p-value distribution from regression analysis
real_plink_output <- read.csv("2020-05-04-108samples-maf005-NOmonomorphic-75support.assoc.fisher",
                         header = TRUE,
                         sep = "")

# adjust real p values
real_plink_output$padj <- p.adjust(real_plink_output$P,
                                   method = "bonferroni")
# vector of names of real sig snps
real_snps <- real_plink_output$SNP[real_plink_output$padj < 0.05]

# compare sig snp numbers


simulation_size <- 1000
sig_snp_num <- c()
sig_snp_names <- list()
real_snps_proportion_in_simulated_sig_snps <- c()

for(ticker in 1:simulation_size){
  # adjust for Bonferroni
  data$file_contents[[ticker]]$Padj <- p.adjust(data$file_contents[[ticker]]$P,
                             method = "bonferroni")
  
  
  adj_pvalue_vec <- data$file_contents[[ticker]]$Padj
  
  # save nb of sig SNPs
  sig_snp_num <- c(sig_snp_num, sum(adj_pvalue_vec < 0.05))
  
  # save name of SNPs
  sig_snp_names[[ticker]] <- data$file_contents[[ticker]]$SNP[data$file_contents[[ticker]]$Padj < 0.05]

  # calculate proportion of real sig snps in simulation
  real_snps_proportion_in_simulated_sig_snps[ticker] <- sum(sig_snp_names[[ticker]] 
                                                       %in% real_snps) / 
    length(sig_snp_names[[ticker]])
  
  }


# make a dataframe
sig_snps_df <- as.data.frame(sig_snp_num)


# make a list of raw p values
raw_p_val_list <- sapply(data$file_contents, `[`, c("P"))
raw_p_val_vec <- unlist(raw_p_val_list)

# make a dataframe
raw_p_val_df <- as.data.frame(raw_p_val_vec)

# plot the distribution of simulated p-values
raw_p_val_df %>%
  ggplot( aes(x=raw_p_val_vec)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha = 0.8)
