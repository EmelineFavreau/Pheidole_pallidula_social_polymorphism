#!/bin/bash
#$ -cwd                        # Run the code from the current directory
#$ -j y                        # Merge the standard output and standard error
#$ -l h_rt=2:00:00             # 1 task - testing if this works
#$ -l h_vmem=4G                # 23G RAM per core (23G total here)
#$ -t 1-2514                    # the number of iterations
#$ -o logs/$JOB_NAME/$JOB_ID/  # redirect the standard output from the job
#$ -tc 100 	                   # take only 100 files at once
INPUT_CONTIG=$(sed -n "${SGE_TASK_ID}p" tmp/5kb_hymenoptera_Ppal_E_contigs)
./calculate-read-depth-for-all-samples.sh $INPUT_CONTIG




#declare -a contig_list= cat input/all_hymenoptera_Ppal_E_contigs


#for contig in "${contig_list[@]}"; do

#  qsub calculate-read-depth-for-all-samples.sh $contig

#done
