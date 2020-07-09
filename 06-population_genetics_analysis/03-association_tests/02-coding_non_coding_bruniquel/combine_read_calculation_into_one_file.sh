#!/bin/bash
contig_list=( contig_1346	contig_1697	contig_1752	contig_1768	contig_2241	contig_237	contig_3784	contig_43	contig_4645	contig_4728	contig_1158	contig_1193	contig_1223	contig_1470	contig_1654	contig_2111)

for contig in "${contig_list[@]}"; do
  paste -s samples > result/${contig}/all-samples-read-count5kb
  paste -d '\t' tmp/${contig}/mean_normalised_by_median* >> result/${contig}/all-samples-read-count5kb

done
