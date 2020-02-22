rm(list=ls())
library(ggplot2)

#get option from command
##############################################################################################
args<-commandArgs(T)
input_file = args[1]
#input_file = "S2_1.forward";input_strand = "-"
#input_file = "S2_1.reverse";input_strand = "+"
if (grepl("forward",input_file)){
  input_strand = "-"
}else{
  input_strand = "+"
}
##############################################################################################

#read the gtf file to get the range for each gene 
gtf = read.csv("../../data/TTS_data/NC_005791_kegg2.gtf",
               sep  = "\t" , stringsAsFactors = F, header = F)
gtf = gtf[gtf$V3=="gene" , ]
gtf$gene = sub(  "gene_id ","" ,gtf$V9)
gtf$gene = sub(";.*","",gtf$gene)
rownames(gtf) <- gtf$gene

#Because strand-specific suquencing was used, it was necessary to distinguish gene in positive and negative strand
pos_strand<- function(x){
  #x = "+"
  tmp = gtf[ gtf$V7==x   ,]
  gtf_for =  tmp
  gtf_for = gtf_for[,c(4,5)]
  gtf_for= apply(gtf_for, 1, function(x) {  c(x[1]:x[2]) } ) 
  gtf_for = stack(gtf_for)
  colnames(gtf_for) <- c( "pos" , "gene")
  return(gtf_for)
}

#
gene_strand<- function(x){
  #x = "-"
  tmp = gtf[ gtf$V7==x   ,]
  gtf_for =  tmp
  gtf_for = gtf_for[,c(4,5)]
  colnames(gtf_for)<-c("start","end")
  gtf_for$gene <-  rownames(gtf_for)
  return(gtf_for)
}

#read the abundance in each locus. 
dat = read.csv("../../data/TTS_data/sample_cov.bed",sep = "\t",stringsAsFactors = F)
#initialize a dataframe in order to store the result
result = read.csv("../../data/TTS_data/TES.csv" ,stringsAsFactors = F )
#select gene in positive or negative strand 
pos_gene = pos_strand(input_strand);pos_gene$gene = as.character(pos_gene$gene)
gene_bed = gene_strand(input_strand)


result_dat = dat[ ,  c( "pos", input_file)]

#add  strand,position,3'end_relative_cov information in result 
result_dat$strand  <- input_strand
result_dat = result_dat[ ,  c("strand" ,"pos" ,input_file  )  ]
colnames(result_dat) <- c("strand","position","X3.end.reads_cov")

#normalize read number(coverage) to relative coverage
result_dat$X3.end_relative_cov = result_dat$X3.end.reads_cov / sum(result_dat$X3.end.reads_cov )

#initialize a dataframe in order to store the result
result_dat$head.reads <-NA
result_dat$X3.end_ratio_.1..1 <-NA
result_dat$X3.end_ratio_.2..2 <-NA

#add location_tag(UTR / ORF) information
result_dat$location_tag.UTR...ORF.  <- ifelse( result_dat$position %in% pos_gene$pos ,
                                               "ORF" , "UTR")

######################add local gene  information
pos_gene <- unstack(pos_gene ,  gene~pos  )
pos_gene <- sapply( pos_gene, 
                    function(x) { paste(x , sep = ";" , collapse = ";") } )

pos_gene = data.frame(gene = unlist(pos_gene) ,pos = names(pos_gene))
pos_gene$gene <- as.character(pos_gene$gene)
rownames( pos_gene ) <- as.character(pos_gene$pos)
rownames(result_dat) <- as.character( result_dat$position)
result_dat$local_gene <- NA
ptm <- proc.time()
result_dat[rownames(pos_gene) , "local_gene"] <-  pos_gene[  , "gene" ]
proc.time() - ptm


######################ad up_gene information
result_dat$upstream.gene <-NA
if ( input_strand == "-" ){
  up_gene = read.csv("../../data/TTS_data/up_gene_neg.csv" , stringsAsFactors = F , header = F)
}else{
  up_gene = read.csv("../../data/TTS_data/up_gene_pos.csv" , stringsAsFactors = F , header = F)
}
rownames(up_gene) <-  up_gene$V3
up_gene$V3<- NULL
up_gene <- apply(up_gene, 1, function(x) { c(x[1]:x[2]) })
up_gene <- stack(up_gene)
colnames(up_gene) = c("pos" , "gene")
#result_dat = result_dat[c(1:20000),];up_gene =  up_gene[  c(1:20000), ]
result_dat$upstream.gene <- up_gene$gene

result_dat$start.codon <-NA
result_dat$stop.codon <-NA
result_dat$dist..to.stop.codon  <-NA
result_dat$sequence.from..50.to..20. <-NA

######################add 3'end_ratio_-1/+1,3'end_ratio_-2/+2 information
ptm <- proc.time()
if ( input_strand == "-" ){
  for(i in 3:( nrow(result_dat)-3  )){
    #i=2
    result_dat[i,"head.reads"] <- result_dat[i,"X3.end.reads_cov"] -result_dat[i-1,"X3.end.reads_cov"]
    result_dat[i,"X3.end_ratio_.1..1"] <- (result_dat[i,"X3.end.reads_cov"] + 0.001) / (result_dat[i-1,"X3.end.reads_cov"] + 0.001)
    result_dat[i,"X3.end_ratio_.2..2"] <- (result_dat[i+1,"X3.end.reads_cov"]+ 0.001) / (result_dat[i-2,"X3.end.reads_cov"] +0.001)
  }
}else{
  for(i in 3:(nrow(result_dat)-3)){
    result_dat[i,"head.reads"] <- result_dat[i,"X3.end.reads_cov"] -result_dat[i+1,"X3.end.reads_cov"]
    result_dat[i,"X3.end_ratio_.1..1"] <- (result_dat[i,"X3.end.reads_cov"] + 0.001)/ (result_dat[i+1,"X3.end.reads_cov"] + 0.001)
    result_dat[i,"X3.end_ratio_.2..2"] <- (result_dat[i-1,"X3.end.reads_cov"]+ 0.001) / (result_dat[i+2,"X3.end.reads_cov"]  + 0.001)
  }
}
proc.time() - ptm


##################add start codon,stop codon information
gene_bed <-  gene_strand(input_strand)
if ( input_strand == "-" ){
  
  result_dat$start.codon = gene_bed[  result_dat$upstream.gene ,"end"  ]
  result_dat$stop.codon = gene_bed[  result_dat$upstream.gene ,"start"  ]
  result_dat$dist..to.stop.codon = result_dat$stop.codon -result_dat$position
}else{
  result_dat$start.codon = gene_bed[  result_dat$upstream.gene ,"start"  ]
  result_dat$stop.codon = gene_bed[  result_dat$upstream.gene ,"end"  ]
  result_dat$dist..to.stop.codon = result_dat$position - result_dat$stop.codon
}

#rounding of numbers
result_dat$X3.end_ratio_.1..1 <- round(result_dat$X3.end_ratio_.1..1,6)
result_dat$X3.end_ratio_.2..2 <- round(result_dat$X3.end_ratio_.2..2,6)
#write output
write.csv( result_dat ,  paste("../../result/TTS_result/",input_file ,".csv" , sep = "") , row.names = F)


