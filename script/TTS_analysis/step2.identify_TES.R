rm(list=ls())
library(dplyr)

#read the gtf file to get the range for each gene 
gtf = read.csv("../../data/TTS_data/NC_005791_kegg2.gtf",
               sep  = "\t" , stringsAsFactors = F, header = F)
gtf = gtf[gtf$V3=="gene" , ]
gtf$gene = sub(  "gene_id ","" ,gtf$V9)
gtf$gene = sub(";.*","",gtf$gene)
rownames(gtf) <- gtf$gene

#get coordinate in 200nt in the downstream of the gene stop codon 
gtf$TES_end = ifelse(gtf$V7=="-" , gtf$V4-200,gtf$V5+200)
gtf$TES_end = ifelse(gtf$TES_end<1,  1,gtf$TES_end)

####sample
args<-commandArgs(T)
input_file = args[1]
####

#input file and initialize some culumn
pos = read.csv(paste0("../../result/TTS_result/",input_file,".reverse3.csv"))
pos [is.na(pos$head.reads),c("head.reads" ,"X3.end_ratio_.1..1" ,"X3.end_ratio_.2..2")] <-c(0,0,0)
pos$position = rownames(pos)
neg = read.csv(paste0("../../result/TTS_result/",input_file,".forward3.csv"))
neg [is.na(neg$head.reads),c("head.reads" ,"X3.end_ratio_.1..1" ,"X3.end_ratio_.2..2")] <-c(0,0,0)
neg$position = rownames(neg)

###############
my_group  <- function(x){
  x = as.numeric(x)
  groupid = 1
  result=c()
  for (i in x){
    if ((i+1) %in% x){
      result=c(result,groupid)
    }else{
      result=c(result,groupid)
      groupid=groupid+1
    }
  }
  return(result)
}


#################################count the TTS  of each gene 
#gtf$TES_pos = NAnrow(gtf)
for (i in 1:nrow(gtf)){
 # i=2
  if (gtf[ i , "V7" ]=="-"){
    TES_start = gtf[ i , "V4" ]
    tmp=neg
  }else{
    TES_start = gtf[ i , "V5" ]
    tmp=pos
  }
  TESrange=c(TES_start: gtf[ i , "TES_end" ])
  TESrange = sort(TESrange)
  tmp = tmp[ TESrange , ]
  tmp = tmp[!is.na(tmp$position),]
  
  tmp1 = tmp[tmp$head.reads > 5 & tmp$X3.end_ratio_.1..1  > 1.1  ,  ]
  tmp2 = tmp[tmp$head.reads > 5 & tmp$X3.end_ratio_.2..2  > 1.1  ,  ]
  
  ##if meet ead.reads > 5(1) and 3.end_ratio_.1..1(2)  > 1.1, the TES_type marked as T1.
  if (nrow(tmp1) > 0){
    tmp = tmp1
    TES_type = "T1"
  } else if(nrow(tmp2) > 0){
    tmp = tmp2
    TES_type = "T2"
  }else{
    gtf[i,"TES_pos_num"] = 0
    gtf[i,"TES_pos"] = NA
    gtf[i,"TES_pos_max"] = NA
    gtf[i,"TES_type"] = "T3"
    next()
  }
  
  #when there were consecutive sites that meet the condition 1and2, the site with highest head.reads was taken as tts
  tmp_ID_group = my_group(tmp$position)
  tmp$ID_group = tmp_ID_group
  tmp_ID_group_max =  by( tmp , tmp$ID_group,function(x){   x$position[which.max(x$head.reads)] } )
  tmp_ID_group_max = as.vector(tmp_ID_group_max)
  
  #add other information
  gtf[i,"TES_pos_num"] = length(tmp_ID_group_max)
  gtf[i,"TES_pos"] = paste(tmp_ID_group_max ,collapse = ";" )
  gtf[i,"TES_pos_max"] = tmp[which.max(tmp$head.reads),"position"]
  gtf[i,"TES_type"] = TES_type
}
##############################################################

##################################add information in - strand 
gtf_neg = gtf[ gtf$V7 == "-" , ]
gtf_neg = gtf_neg[ ! is.na(gtf_neg$TES_pos) , ]
for ( i in 1:nrow(gtf_neg)){
 # i=141
  pos1= gtf_neg[i,"TES_pos"] %>% strsplit(split = ";") %>% unlist
  pos2= gtf_neg[i+1,"TES_pos"] %>% strsplit(split = ";") %>% unlist
  if (length( intersect(pos1 ,pos2) ) !=0){
    pos1 = pos1[!pos1 %in% pos2]

   gtf_neg[i,"TES_pos_num"] = length(pos1)
   TES_pos_max = neg[ c(pos1) , ]
   TES_pos_max = TES_pos_max$position[which.max(TES_pos_max$head.reads)]
  
  if  (length(pos1) == 0){
    gtf_neg[i,"TES_pos_max"] = NA
    gtf_neg[i,"TES_pos"] = NA
  }else{
    gtf_neg[i,"TES_pos_max"] = TES_pos_max
    pos1 = paste(pos1 , collapse = ";")
    gtf_neg[i,"TES_pos"] = pos1
  }
  }
  
} 
tmp_gtf_neg = gtf_neg
##############################################################

##################################add information in + strand 
gtf_neg = gtf[ gtf$V7 == "+" , ]
gtf_neg = gtf_neg[ ! is.na(gtf_neg$TES_pos) , ]
for ( i in 1:nrow(gtf_neg)){
  # i=141
  pos1= gtf_neg[i,"TES_pos"] %>% strsplit(split = ";") %>% unlist
  pos2= gtf_neg[i+1,"TES_pos"] %>% strsplit(split = ";") %>% unlist
  if (length( intersect(pos1 ,pos2) ) !=0){
    pos2 = pos2[!pos2 %in% pos1]
    
    gtf_neg[i+1,"TES_pos_num"] = length(pos2)
    TES_pos_max = neg[ c(pos2) , ]
    TES_pos_max = TES_pos_max$position[which.max(TES_pos_max$head.reads)]
    
    if  (length(pos2) == 0){
      gtf_neg[i+1,"TES_pos_max"] = NA
      gtf_neg[i+1,"TES_pos"] = NA
    }else{
      gtf_neg[i+1,"TES_pos_max"] = TES_pos_max
      pos2 = paste(pos2 , collapse = ";")
      gtf_neg[i+1,"TES_pos"] = pos2
    }
  }
  
} 
tmp_gtf_pos = gtf_neg
##############################################################

##################################combine result and output 
tmp_gtf = rbind(tmp_gtf_neg , tmp_gtf_pos , gtf[ is.na(gtf$TES_pos), ])
tmp_gtf = tmp_gtf[ order(tmp_gtf$V4) , ]
write.csv(x = tmp_gtf ,
          paste0("../../result/TTS_result/",input_file,"_TES.gtf.csv_tmp"),  
          row.names = F)

#############################recheck TES
#rm(list=ls())

dat = read.csv( paste0("../../result/TTS_result/",input_file,"_TES.gtf.csv_tmp"),
               stringsAsFactors = F)
rownames(dat) <- dat$gene
dat  = dat[ !dat$TES_type == "T3" , ]
dat  = dat[dat$TES_pos_num==0 , ]


pos = read.csv(paste0("../../result/TTS_result/",input_file,".reverse3.csv")) 
pos [is.na(pos$head.reads),c("head.reads" ,"X3.end_ratio_.1..1" ,"X3.end_ratio_.2..2")] <-c(0,0,0)
pos$position = rownames(pos)
neg = read.csv(paste0("../../result/TTS_result/",input_file,".forward3.csv"))
neg [is.na(neg$head.reads),c("head.reads" ,"X3.end_ratio_.1..1" ,"X3.end_ratio_.2..2")] <-c(0,0,0)
neg$position = rownames(neg)
my_group  <- function(x){
  x = as.numeric(x)
  groupid = 1
  result=c()
  for (i in x){
    if ((i+1) %in% x){
      result=c(result,groupid)
    }else{
      result=c(result,groupid)
      groupid=groupid+1
    }
  }
  return(result)
}

gtf = dat

#gtf$TES_pos = NAnrow(gtf)
for (i in 1:nrow(gtf)){
  # i=1 
  if (gtf[ i , "TES_type" ]=="T2"){
    gtf[ i , "TES_type" ]="T3"
    next()
  }
  
  if (gtf[ i , "V7" ]=="-"){
    TES_start = gtf[ i , "V4" ]
    tmp=neg
  }else{
    TES_start = gtf[ i , "V5" ]
    tmp=pos
  }
  TESrange=c(TES_start: gtf[ i , "TES_end" ])
  TESrange = sort(TESrange)
  tmp = tmp[ TESrange , ]
  tmp = tmp[!is.na(tmp$position),]
  
  
  tmp2 = tmp[tmp$head.reads > 5 & tmp$X3.end_ratio_.2..2  > 1.1  ,  ]
  tmp2 = tmp2[tmp2$head.reads > 5 & tmp2$X3.end_ratio_.1..1  < 1.1  ,  ]
  
  if(nrow(tmp2) > 0){
    tmp = tmp2
    TES_type = "T2"
  }else{
    gtf[i,"TES_pos_num"] = 0
    gtf[i,"TES_pos"] = NA
    gtf[i,"TES_pos_max"] = NA
    gtf[i,"TES_type"] = "T3"
    next()
  }
  
  
  tmp_ID_group = my_group(tmp$position)
  tmp$ID_group = tmp_ID_group
  tmp_ID_group_max =  by( tmp , tmp$ID_group,function(x){   x$position[which.max(x$head.reads)] } )
  tmp_ID_group_max = as.vector(tmp_ID_group_max)
  gtf[i,"TES_pos_num"] = length(tmp_ID_group_max)
  gtf[i,"TES_pos"] = paste(tmp_ID_group_max ,collapse = ";" )
  gtf[i,"TES_pos_max"] = tmp[which.max(tmp$head.reads),"position"]
  gtf[i,"TES_type"] = TES_type
}

dat = read.csv(paste0("../../result/TTS_result/",input_file,"_TES.gtf.csv_tmp"), 
               stringsAsFactors = F)
rownames(dat) <- dat$gene
dat[ rownames(gtf) , ] <- gtf

write.csv(x = dat, 
          paste0("../../result/TTS_result/",input_file,"_TES.gtf.csv"), 
          row.names = F)
file.remove(paste0("../../result/TTS_result/",input_file,"_TES.gtf.csv_tmp"))
