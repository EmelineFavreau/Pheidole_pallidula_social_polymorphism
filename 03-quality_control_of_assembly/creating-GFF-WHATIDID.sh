module load  augustus/3.2.3
module load  maker/2.31.9
maker -CTL

#edited; remove repeatmasking because it was crashing.                                                                                                                                          
# made it only use solenopsis prots as hints - nothing else (!) and nasonia augustus. just to get sth crap v v quick. 
 
# qlogined to sm10; ran many small jobs (of 4 cores) because maker spins jobs up and down so rapidly that it wasn't using much of the 48 cores available. 

maker -cpu 4
maker -cpu 4
maker -cpu 4
maker -cpu 4
maker -cpu 4
maker -cpu 4
maker -cpu 4
maker -cpu 4
maker -cpu 4
maker -cpu 4
maker -cpu 4
maker -cpu 4
maker -cpu 4
maker -cpu 4
maker -cpu 4
maker -cpu 4

find .  -type f -name '*gff' ! -path '*theVoid*' | parallel cat {} > reference.gff
tar -czvf reference.maker.output.tgz reference.maker.output         tar -czvf reference.maker.output.tgz reference.maker.output         
tar -czvf reference.maker.output.tgz reference.maker.output
