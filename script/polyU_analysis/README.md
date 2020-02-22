# This script was used to count poly(U) around the gene transcription termination site.

## The input files required for this script：
        species.fasta
        species.gff3

## The software required for this script：      
        samtools
  
### Each run:
```
    Rscript  polyU.R  -f  "../../data/polyU_data/Acidilobus saccharovorans/Acidilobus saccharovorans genome.fasta"  -g "../../data/polyU_data/Acidilobus saccharovorans/sequence.gff3"  -o  '../../result/polyU_result/Acidilobus saccharovorans'  
```

### Batch run:
```
   Rscript batch_polyU.R
```