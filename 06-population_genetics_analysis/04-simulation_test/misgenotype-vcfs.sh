#!/bin/bash

# Create 1000 VCFs with 5% misgenotyped samples
# This is part of a test of analysis power

# Copyright Emeline Favreau, Queen Mary University of London

### Objectives:
# Take the original dataset VCF
# Assign different samples to the real VCF for 5%
# Re-run Fisher

### Project preparation ###
module load bcftools/1.8
module load vcftools/0.1.15
module load plink/1.9-170906



for ticker in {1..1000}
do



# 5% samples are assigned the other social type
# misgenotyping model with Fisher's exact test and Bonferroni adjustment
todayanalysis="2021-01-13-misgenotyping-5-Fisher-Bonferroni${ticker}"

# Step 1: Change the social type of 5% of samples to the VCF file

# Make misgenotyping in sample file

# number of samples: 108
# wc -l result/sample_names.txt

# number that I should change: 5 (= 5%)

# sample randomly that number
shuf -n 5 result/sample_names.txt > tmp/${todayanalysis}_sample_names_to_misgenotype.txt

# change social type in a temp file
grep "P" tmp/${todayanalysis}_sample_names_to_misgenotype.txt \
	| sed 's/-P/-M/g' > tmp/${todayanalysis}_sample_names_misgenotyped.txt

grep "M" tmp/${todayanalysis}_sample_names_to_misgenotype.txt \
	| sed 's/-M/-P/g' >> tmp/${todayanalysis}_sample_names_misgenotyped.txt

grep "N" tmp/${todayanalysis}_sample_names_to_misgenotype.txt \
	| sed 's/-N/-P/g' >> tmp/${todayanalysis}_sample_names_misgenotyped.txt

# sort resulting files
# original
sort tmp/${todayanalysis}_sample_names_to_misgenotype.txt \
	> tmp/${todayanalysis}_original

# altered
sort tmp/${todayanalysis}_sample_names_misgenotyped.txt \
	> tmp/${todayanalysis}_altered

# combine into one file that takes two columns
paste tmp/${todayanalysis}_original tmp/${todayanalysis}_altered \
	> tmp/${todayanalysis}_patterns.txt



# Step 2: make misgenotyping in VCF (121,000 SNPs)

# make a new file
cp input/2020-05-04-108samples-maf005-NOmonomorphic-75support.recode.vcf \
	tmp/${todayanalysis}.vcf

# replace in line in the VCF
while read -r pattern replacement; do   
    sed -i "s/$pattern/$replacement/" tmp/${todayanalysis}.vcf
done < tmp/${todayanalysis}_patterns.txt



# Step 3: change pheno.txt
cp pheno.txt tmp/${todayanalysis}_pheno.txt

while read -r pattern replacement; do   
    sed -i "s/$pattern/$replacement/g" tmp/${todayanalysis}_pheno.txt
done < tmp/${todayanalysis}_patterns.txt

# change social type code in phenotype column
echo -e "#FID\tIID\tphenotype" > tmp/${todayanalysis}_pheno_cleaned.txt

grep "P" tmp/${todayanalysis}_pheno.txt \
	| sed 's/2$/1/g' >> tmp/${todayanalysis}_pheno_cleaned.txt


grep "\-M" tmp/${todayanalysis}_pheno.txt \
	| sed 's/1$/2/g' >> tmp/${todayanalysis}_pheno_cleaned.txt


grep "\-N" tmp/${todayanalysis}_pheno.txt \
	| sed 's/2$/1/g' >> tmp/${todayanalysis}_pheno_cleaned.txt



# Step 4: combine real SNPs with the simulated SNPs
plink --vcf tmp/${todayanalysis}.vcf \
       --allow-extra-chr \
       --allow-no-sex \
       --pheno tmp/${todayanalysis}_pheno_cleaned.txt \
       --set-missing-var-ids @:#\$1,\$2 \
       --make-bed \
       --out tmp/${todayanalysis}-id2

# recode as VCF file
plink --bfile tmp/${todayanalysis}-id2 \
       --allow-extra-chr \
       --allow-no-sex \
       --recode vcf-fid  \
       --out tmp/${todayanalysis}-updated

# obtain SNP matrix for R script of genotype likelihood test
bcftools query -f '%ID\t%POS[\t%GT]\n' tmp/${todayanalysis}-updated.vcf \
      > tmp/${todayanalysis}-snp_matrix.txt

# obtain sample vec
bcftools query -l tmp/${todayanalysis}-updated.vcf \
             > tmp/${todayanalysis}-sample_names.txt




# Step 5: run a Fisher test of association for social type
# standard case/control association analysis using Fisher's exact test to generate significance
plink --vcf tmp/${todayanalysis}-updated.vcf \
      --allow-extra-chr \
      --allow-no-sex \
      --pheno tmp/${todayanalysis}_pheno_cleaned.txt \
      --assoc fisher \
      --out tmp/${todayanalysis}

done