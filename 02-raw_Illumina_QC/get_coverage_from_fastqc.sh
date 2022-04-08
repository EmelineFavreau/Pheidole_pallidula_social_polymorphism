#!/bin/bash
while read -r LINE; do
	unzip ../2018-05-30-read-qc/results/"$LINE"_fastqc.zip
	echo "$LINE" >> 115_pheidole_coverage.txt
	cd "$LINE"_fastqc
	sed '7q;d' fastqc_data.txt >> ../115_pheidole_coverage.txt
	cd ../
done < list_of_samples_0001.txt



