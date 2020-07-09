#!/bin/bash
# extract uniref IDs

INPUTFILE=$1
OUTPUTFILE=$2

## get uniref ids and put into urls
cut -f 2 $INPUTFILE | sort | uniq | parallel echo https://www.uniprot.org/uniref/{}.fasta > urls.txt
	
## download fasta from urls
cat urls.txt | parallel wget -q -O - | grep ">" > $OUTPUTFILE

## clean up
rm urls.txt
