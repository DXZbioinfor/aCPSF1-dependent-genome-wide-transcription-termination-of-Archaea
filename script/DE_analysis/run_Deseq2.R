
#### The input file named readcount.csv was provided in supplementary materials. The output file named DF_gene.csv was used in follow analysis.

library("DESeq2")
rm(list=ls())

# input pheno data 
pheno_data = data.frame(Diagonsis = c(rep("S2",3) , rep("CPSF1",3) )  )
pheno_data$Diagonsis = as.factor(pheno_data$Diagonsis)

# input readcount data  
mRNA_expre = read.csv("../../data/TES_data/readcount.csv" ,header = T,row.names = 1)
KEGGID<-mRNA_expre$KEGGID
mRNA_expre$KEGGID  <- NULL
rownames(pheno_data) = colnames(mRNA_expre)

# transfrom to Deseq2Matrinx
dds_mRNA = DESeqDataSetFromMatrix(countData = mRNA_expre  ,
                                    colData = pheno_data,  
                                    design =~ Diagonsis ) 

# estimate mRNA sizefactor
dds_mRNA = estimateSizeFactors(dds_mRNA)   

# calculate different expression gene 
dds2 <- DESeq(dds_mRNA,parallel=TRUE,betaPrior=FALSE)
res <-  results(dds2, contrast=c("Diagonsis","CPSF1","S2"),parallel=TRUE)

# output result 
res_dataframe=as.data.frame(res)
res_dataframe=data.frame(ID = rownames(res_dataframe),GeneID=KEGGID,res_dataframe)
res_dataframe=res_dataframe[order(res_dataframe$GeneID),]
write.csv(res_dataframe,"../../result/DE_result/DE_gene.csv",row.names = F)

# report 
res_dataframe = res_dataframe[ !is.na(res_dataframe$padj) , ]
cat("Number of differentially expressed genes(adjp<0.05): ", table(res_dataframe$padj< 0.05)[2])

cat("Number of Up-regulated genes: ",nrow(res_dataframe[res_dataframe$log2FoldChange>0 & res_dataframe$padj < 0.05,] ))
cat("Number of Down-regulated genes: ",nrow(res_dataframe[res_dataframe$log2FoldChange<0 & res_dataframe$padj < 0.05,] ))

