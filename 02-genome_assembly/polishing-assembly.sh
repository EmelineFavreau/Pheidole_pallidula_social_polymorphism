# using Pilon to polish # uses template.sh
# All computations were done on srm nodes.


# First address SNP and small indels.
for i in $(seq 1 10); do
  mkdir pilon
  params="$(( $i - 1 )).asm.fa R1.fq.gz R2.fq.gz snps,indels" ./polish.sh
  mv pilon pilon.${i}; ln -s pilon.${i}/pilon.fasta ${i}.asm.fa
done

# Try ambiguous bases and gap-filling.
for i in $(seq 11 20); do
  mkdir pilon
  params="$(( $i - 1 )).asm.fa R1.fq.gz R2.fq.gz snps,indels,amb,gaps" ./polish.sh
  mv pilon pilon.${i}; ln -s pilon.${i}/pilon.fasta ${i}.asm.fa
done

