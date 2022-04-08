# assessing sequence length and nucleotide composition of PASS basecalled MinION reads
module load seqtk

# get average length # Output format: chr, length, #A, #C, #G, #T, #2, #3, #4,
#CpG, #tv, #ts, #CpG-ts
cat input/ppal16/*.fastq | seqtk comp  > result/ppal16-composition
cat input/ppal15/*.fastq | seqtk comp  > result/ppal15-composition
cat input/ppal4/*.fastq | seqtk comp  > result/ppal4-composition
cat input/phei2/*.fq | seqtk comp  > result/phei2-composition
cat input/phei3/*.fq | seqtk comp  > result/phei3-composition
cat input/phei4/*.fq | seqtk comp  > result/phei4-composition

# quick calculations for each run
minion_composition.Rmd
