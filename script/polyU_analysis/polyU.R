rm(list = ls())
library(optparse)
library(dplyr)


#get options from command 
option_list <- list( 
    make_option(c("-f", "--fasta"), type="character",help="Peference genome for sepcies"),

    make_option(c("-g", "--gff3"), type="character",help="Gff3 file for sepcies"),
    
    make_option(c("-o", "--outfile"), type="character",help="Prefix for output file" )

    )
opt <- parse_args(OptionParser(option_list=option_list))

options(scipen=200)

#read gff file and  selecct the line for  gene 
dat = read.csv(opt$gff3, header = F , sep="\t",comment.char = "#")
dat = dat[dat$V3=="gene",]
BED=dat

#get the sequence around the transcription termination site
my_seq <- function( chr, s, pos_start, pos_end){
    fasta = opt$fasta
    if (s=="+"){
        conmand=sprintf("~/bin/samtools   faidx  \"%s\" %s:%s-%s",fasta,chr ,pos_start , pos_end  )
    } else{
        conmand=sprintf("~/bin/samtools   faidx   -i  \"%s\" %s:%s-%s",fasta,chr ,pos_start , pos_end )
    }
    a  <- system(conmand,intern = T)
    return( a )
    
}

#count the polyT in sequence 
my_polyT <- function(x){
    #x<- result[i,"up_seq"]
    x <- gsub("A|C|G" , " " , x)
    x <- strsplit(x, split = " ")[[1]]
    x <- x[x!=""]
    count_T <- table(x)
    
    #if ployT>7 classfied  to  ployT7
    T_times= c(8:30)
    count_T_greater7 <- lapply(T_times , function(x) {  paste(rep("T" , x) ,collapse = "" ) } )
    count_T_greater7 = lapply(count_T_greater7, function(x){ count_T[x] } )
    count_T_greater7 = count_T_greater7 %>% unlist %>% na.omit %>% sum
    
    #cout for ployT4-7
    count_T <- c(count_T["TTTT"] ,count_T["TTTTT"]  , count_T["TTTTTT"]  , count_T["TTTTTTT"] +count_T_greater7   )

    count_T[is.na (count_T)] <- 0
    a = paste(count_T , collapse  = ";")
    b= sum(count_T)
    return(c(a,b))
}

#initialize a dataframe in order to store the result
result = data.frame( BED$V1 , BED$V4 , BED$V5, BED$V7 , BED$V9 )
colnames(result)= c("chr" , "start"  , "end" , "strand" , "gene")
result$up_name = NA;result$up_seq = NA ; result$up_T4567 <-  NULL  ;result$up_T4567_total <-  NULL 
result$down_name = NA ;result$down_seq = NA ; result$down_T4567 <-  NULL  ;result$down_T4567_total <-  NULL 

#add information for each gene
for (i in c(1:nrow(dat))){
    #i=411
    chr=BED[i,1] %>% as.character()
    s=BED[i , 7]  %>% as.character()
    if ( s == "+"){ 
        up_start = BED[i,5]-199 ; up_end = BED[i,5]
        stram_start= BED[i,5] ;  stram_end = BED[i,5] + 199 
        up_start=ifelse(up_start<0 , 1, up_start)
    }else{
        up_start = BED[i,4] ;  up_end= BED[i,4] + 199   
        stram_end=  BED[i,4] ; stram_start = BED[i,4] - 199
        stram_start=ifelse(stram_start<0 , 1, stram_start)
    }
    up_seq=my_seq(chr = chr,s = s,up_start,  up_end)
    stram_seq = my_seq(chr = chr,s = s,stram_start,  stram_end)
    result[i,"up_name"] <- up_seq[1]
    result[i,"up_seq"] <- paste(up_seq[2:5] %>% na.omit() ,collapse = "")
    result[i,"up_T4567"] <- my_polyT(    result[i,"up_seq"] )[1]
    result[i,"up_T4567_total"] <- my_polyT(    result[i,"up_seq"] )[2]
    
    result[i,"down_name"] <- stram_seq[1]
    result[i,"down_seq"] <- paste(stram_seq[2:5]  %>% na.omit()  ,collapse = "")   
    result[i,"down_T4567"] <- my_polyT(    result[i,"down_seq"] )[1]
    result[i,"down_T4567_total"] <- my_polyT(    result[i,"down_seq"] )[2]
    
    
}

#add species in each line
result$species=opt$outfile
write.csv(result , file = paste0(opt$outfile , ".polyU.csv") , row.names = F)
