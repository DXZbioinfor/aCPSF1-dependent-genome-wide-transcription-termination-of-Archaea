# Those codes were used for identify transcription termination sites and draw metapictures

## Step1. Scan each locus
We scanned each locus in genome to get the basic information in each locus such as abundance，coverage，upstream and downstream gene.etc. S2_1 and DCL_2 represent wild type and mutant type, respectively. Forward and reverse represent - and + strand, respectively. 

Run by:
``` 
Rscript    step1.all_site.R  S2_1.forward
Rscript    step1.all_site.R  S2_1.reverse	
Rscript    step1.all_site.R  DCL_2.forward	
Rscript    step1.all_site.R  DCL_2.reverse
```

The output files were used as input file in step2
```
DCL_2.forward3.csv
DCL_2.reverse3.csv
S2_1.forward3.csv
S2_1.reverse3.csv
```


## step2. Detect candidate TTSs	
We used four criteria to detect TTSs. (1) within 200 nt of the downstream region distant to the stop codon of a gene; (2) >1.1 reads ratio of -1 site (predicted TTS) to +1 site (downstream TTS); (3) read-counts of -1 site minus +1 site > 5; (4) upon (2-3) satisfied, the site that has the highest read counts of -1 site minus +1 site and occurs in each of the two biological replicates was recorded as a primary TTS.

Run by:
```
Rscript    step2.identify_TES.R  DCL_2
Rscript    step2.identify_TES.R  S2_1
```

The output files are:
```
S2_1_TES.gtf.csv
DCL_2_TES.gtf.csv
```
In output file, each line records the information about all TSS of a gene. The genes that passed the criteria(1,2,3) were marked as T1 in the column of TES_type. The coordinate marked in the column of *TES_pos_max* was the highest "head reads" locus and was the primary TTS.

## Step3. Convert format
The  *S2_1_TES.gtf.csv* was converted format and add more information.

Run by:
```
Rscript     step3.moreInfor_TES.R
```

The  output file is:
```
TES_information.csv
```


Step4. Plot
Plot the abundance in 20bp around primary TTS in wild and mutant type

Run by:
```
Rscript     step4_metapicture.R
```

The output file is:
```
 metapicture.pdf
```
![image](https://github.com/DXZbioinfor/aCPSF1-dependent-genome-wide-transcription-termination-of-Archaea/blob/master/result/TTS_result/metapicture.png)
