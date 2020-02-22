rm(list=ls())
library("dplyr")

#input data and filter gene with no TTS
dat = read.csv(  paste0("../../result/TTS_result/","S2_1","_TES.gtf.csv"), 
               stringsAsFactors = F)
rownames(dat) <- dat$gene
dat  = dat[ !dat$TES_type == "T3" , ]

#spilte TTS by gene
dat_list = split(dat[ , c("TES_pos" )], dat$gene)
dat_list = sapply(dat_list, function(x){
                        strsplit(x[[1]],split = ";")}
                  )
dat_dat = stack(dat_list)
colnames(dat_dat) <- c("pos","gene")
dat_dat$strand = dat[  dat_dat$gene %>% as.character(), "V7"]
dat_dat$type = dat[  dat_dat$gene %>% as.character(), "TES_type"]

################################################Supplement the basic information of each TTS in the gene according to to the strand where the gene is located
my_peak <- function(x){
 # x = "-"
  #########################main  TES###########################
  peak = dat_dat[ dat_dat$strand==x , ]
  tmp_dat <- dat[dat$V7 == x ,]
  #yes represent there is the primary TTS
  peak$peak <- ifelse(peak$pos %in% tmp_dat$TES_pos_max,"yes","no")
  
  #########################other information###########################
  if (x=="-"){ 
      S2 = read.csv("../../result/TTS_result/S2_1.forward3.csv")
      DCL = read.csv("../../result/TTS_result/DCL_2.forward3.csv")
  }else{
      S2 = read.csv("../../result/TTS_result/S2_1.reverse3.csv")
      DCL = read.csv("../../result/TTS_result/DCL_2.reverse3.csv")    
  }
  base_info = S2[ peak$pos %>% as.character() %>% as.numeric(),  ]
  peak = cbind(peak ,
               base_info[ , 
                          c("location_tag.UTR...ORF.", "local_gene","upstream.gene",
                            "start.codon" ,"stop.codon", "dist..to.stop.codon" )
                        ]
               )
  S2_abun =  S2[ peak$pos %>% as.character() %>% as.numeric(),  c("head.reads", "X3.end_ratio_.1..1"  ,  "X3.end_ratio_.2..2" )]
  colnames(S2_abun) <- paste("S2",c("head.reads"  ,  "end_ratio11", "end_ratio22"),sep = "_")
  
  DCL_abun =  DCL[ peak$pos %>% as.character() %>% as.numeric(),  c("head.reads", "X3.end_ratio_.1..1"  ,  "X3.end_ratio_.2..2" )]
  colnames(DCL_abun) <- paste("DCL",c("head.reads"  ,  "end_ratio11", "end_ratio22"),sep = "_")
  
  peak = cbind(peak,S2_abun, DCL_abun )
  
  #########################result ###########################
  return(peak)
}
################################################
dat_dat = rbind( my_peak("+") , my_peak("-") )


################################################ get the range that need estract seuqnece 
Tes_range <-  function(x,pos,up,down){
  #x= 100;up=20;down=20
  pos = as.numeric(pos)
  if (x=="-"){ 
      Tes_up=pos+up
      Tes_down=pos-down
  }else{
      Tes_up=pos-up
      Tes_down=pos+down  
  }
  return(sort(c(Tes_up ,Tes_down)))
}

################################################ get sequence by samtools 
my_sequence <- function(x,pos=pos,up=up,down=down){
    myrange <- Tes_range (x,pos=pos,up=up,down=down )
    if(x=="-"){
      comand = sprintf("samtools   faidx  -i  ../../data/TTS_data/genome.fasta  BX950229.1:%s-%s",myrange[1],myrange[2])
    }else{
      comand = sprintf("samtools   faidx     ../../data/TTS_data/genome.fasta BX950229.1:%s-%s",myrange[1],myrange[2])
    }
    sequence  <- system(comand,intern = T)
    sequence <- sequence[-1]
    sequence <- paste(sequence,collapse = "")
  
}


dat_range <- dat_dat[,c(1,3)]

for (i in 1:nrow(dat_range)){
  #i=2
  x=dat_range[i ,"strand"];pos=dat_range[i ,"pos"]
  dat_dat$sequence_50_20[i] <-  my_sequence( x,pos=pos,up=50,down=20)
  dat_dat$sequence_30_20[i] <-  my_sequence( x,pos=pos,up=30,down=20)  
  dat_dat$sequence_50_50[i] <-  my_sequence( x,pos=pos,up=50,down=50)  
}
################################################

dat_dat = dat_dat[ order(dat_dat$pos %>% as.numeric()) , ]
write.csv(x = dat_dat , 
          file = "../../result/TTS_result/TES_information.csv",
          row.names = F,quote = F)
