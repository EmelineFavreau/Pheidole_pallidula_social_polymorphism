# obtain N50 from various Hymenoptera 
esearch -db assembly -query 'txid7399[organism:exp]' \
  | esummary \
  | xtract -pattern DocumentSummary -element AssemblyAccession,AssemblyName,SpeciesTaxid,Organism,ContigN50,Scaffold  N50 > hymenoptera-N50-date-stats

