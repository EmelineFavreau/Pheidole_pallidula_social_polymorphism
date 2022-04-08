# Basecalling six MinION runs

# Library from pooled Phei2, Phei3, Phei4:
#basecalling in the cloud (Metrichor)

# Library Ppal4
. ~/bin/activate
{time ~/bin/read_fast5_basecaller.py -r --flowcell FLO-MIN106 --kit SQK-LSK108 --output_format fast5,fastq --input Ppal4/ --save_path output --worker_threads 10 ; } 2> ppal4-time.txt
deactivate

# Library Ppal15
. ~/bin/activate
{time ~/bin/read_fast5_basecaller.py -r --flowcell FLO-MIN106 --kit SQK-LSK108 --output_format fastq --input ppal15/ --save_path ppal15_output --worker_threads 20 ; } 2> ppal15-time.txt
deactivate

# Library Ppal16
. ~/bin/activate
{time ~/bin/read_fast5_basecaller.py -r --flowcell FLO-MIN106 --kit SQK-LSK108 --output_format fastq --input ppal16/ --save_path ppal16_output --worker_threads 5 ; } 2> ppal16-time.txt
deactivate

