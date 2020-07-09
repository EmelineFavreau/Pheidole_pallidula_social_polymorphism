#!/bin/bash
mkdir -p result/read_coverage_per_5kb_all_samples

for contig in $(cat tmp/5kb_hymenoptera_Ppal_E_contigs)
  do
    paste -s samples-without-E15 > result/read_coverage_per_5kb_all_samples/${contig}_all-samples-read-count5kb
    paste -d '\t' tmp/${contig}/mean_normalised_by_median* >> result/read_coverage_per_5kb_all_samples/${contig}_all-samples-read-count5kb

    done
