                                                                         
module load python/2.7.15
mkdir -p tmp/agouti
python ../../soft/AGOUTI/agouti.py scaffold -assembly input/reference.fa -bam input/rnaseq_bwa.bam -gff input/reference.gff -outdir tmp/agouti > agouti.log



# stats are in log file
mv tmp/agouti/agouti* results 


#This improved N50 by 30% - see stats below

#Some of the raw data is in:  
#    results/agouti_denoise/agouti.agouti.join_pairs.noise_free.txt
#However, it does some additional magic to eliminate the unlikely ones, such as Ppal_E.contig_1349 with Ppal_E.contig_1365 (which would imply that the gene spans ~1.5 megabase), and Ppal_E.contig_875 linked to Ppal_E.contig_882.  

#The overall effects (orientations, etc, as would be useful for something like allmaps) are in 
#   results/agouti.agouti.scaffolding_paths.txt
#A different format (which includes amount of support - useful for filtering) is in: 
#   results/agouti.agouti_scaffolding.graph.dot


#Output
#- number of contigs scaffoled: 318
#- number of scaffolds: 142
#- number of contigs in the final assembly: 3954
#- Final assembly N50: 587,760

# input was: 4130 contigs, with assembly N50: 446424
