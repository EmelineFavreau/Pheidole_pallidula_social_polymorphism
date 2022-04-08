#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=4:00:00
#$ -t 1-230

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" list_of_files.txt)
module load fastqc/0.11.5
fastqc -o /data1/autoScratch/monthly/btx077/2018-05-115-pheidole-gwas/results/2018-05-30-read-qc/results/ --noextract /data1/autoScratch/monthly/btx077/2018-05-115-pheidole-gwas/results/2018-05-fastqc-raw-data/input/$INPUT_FILE

