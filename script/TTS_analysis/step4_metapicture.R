rm(list=ls())
library(rlist)
library(ggplot2)
library(reshape2)
library(dplyr)

#input S2 TTS dataframe and select the primary TTS`
dat = read.csv("../../result/TTS_result/TES_information.csv",stringsAsFactors = F)
dat = dat[ dat$peak == "yes"  , ]

#Determine the coordinates to be plot
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


#Extract the abundance information for each site that around primary TTS 20nt in S2 and aCPSF1 sample
my_peak_range_abundance <- function(x){
  # x = "+"
  #########################main  TES###########################
  peak = dat[ dat$strand==x , ]
  #########################headreads  and  relative abundance###########################
  if (x=="-"){ 
    S2 = read.csv("../../result/TTS_result/S2_1.forward3.csv")
    DCL = read.csv("../../result/TTS_result/DCL_2.forward3.csv")
  }else{
    S2 = read.csv("../../result/TTS_result/S2_1.reverse3.csv")
    DCL = read.csv("../../result/TTS_result/DCL_2.reverse3.csv")       
  }
  #########################abundance of range  20-20 ###########################
  result_strain = list()
  for (i in 1:nrow(peak)){
        #i=1
        pos<-peak[i,"pos"]
        my_range<-Tes_range(x,pos = pos,up = 20 ,down = 20 )
        
        S2_abu = S2[ my_range[1]:my_range[2],c("X3.end_relative_cov" , "head.reads","X3.end_ratio_.1..1","X3.end.reads_cov")]
        S2_abu$gene = peak[i,"gene"]
        S2_abu$type = "S2"
      
        DCL_abu = DCL[ my_range[1]:my_range[2],c("X3.end_relative_cov" , "head.reads","X3.end_ratio_.1..1","X3.end.reads_cov")]
        DCL_abu$gene = peak[i,"gene"]
        DCL_abu$type = "DCL"  
        
        if (x=="-"){ 
            S2_abu$ID = c(41:1)
            DCL_abu$ID = c(41:1)
        }else{
            S2_abu$ID = c(1:41)
            DCL_abu$ID = c(1:41) 
        }  
        S2_abu$range_abundance_relative <-  S2_abu$X3.end.reads_cov/sum(S2_abu$X3.end.reads_cov)
        DCL_abu$range_abundance_relative <-  DCL_abu$X3.end.reads_cov/sum(DCL_abu$X3.end.reads_cov)
        
        abu = rbind(S2_abu , DCL_abu)
        abu$strain = x
        result_strain[[peak[i,"gene"]]] = abu
    }
  #########################result ###########################
  return( result_strain)
}

a = my_peak_range_abundance("+") 
b =  my_peak_range_abundance("-") 
dat_dat = rbind(list.rbind(a ),list.rbind(b) ) %>% as.data.frame()


##################################################polt  
#normalize reads of 20 nt each upstream and downstream of all the primary TTSs that were normalized via dividing the individual average read of each site from -21 to +20 by that of -21 site 
dat <- dat_dat
my_norm_abun <- function(x){
 # x = "DCL"
  tmp2 = dat[ dat$type == x , ]
  tmp2 =split(tmp2 , f = tmp2$gene)
  for (i in tmp2){
     gene = i
     start_ab =  gene[,"X3.end.reads_cov"][gene$ID==1]
     gene$X3.end.reads_cov_percent =  gene$X3.end.reads_cov/start_ab
     tmp2[[gene[1,"gene"]]]=gene
  }
  tmp2 = list.rbind(tmp2) %>% as.data.frame() 
  tmp = tmp2  
  ###############################################
  tmp = group_by(tmp ,ID )
  #calculate the average normalized abundance of those sites which have same position relative to respective primary TTS
  tmp = summarise( tmp, 
                   my_percent = mean(X3.end.reads_cov_percent,na.rm = T),
                   my_percent_sd = sd(X3.end.reads_cov_percent,na.rm = T),
                   my_range_abundance_relative =  mean(range_abundance_relative,na.rm = T),
                   my_sd_abundance = sd(range_abundance_relative,na.rm = T))
  tmp =  as.data.frame(tmp)
  tmp$type = x
  return(tmp)

}
dat = rbind(my_norm_abun("S2"), my_norm_abun("DCL") ) 

#melt data and plot picture 
dat<-melt(dat, measure.vars = c("my_percent" ))
p<-ggplot(dat,aes(x=ID,y =value*100 ))+
  geom_step(aes(color = type),stat = 'identity')+ 
  theme_classic(12)+
  theme(aspect.ratio = 1/1.6,legend.title = element_blank(),
        legend.position = c(0.9,0.9),axis.title = element_text(size=10))+
  scale_x_continuous(breaks = seq(1,41))+
  scale_color_manual(values = c("#b2182b","#023858"))+
  scale_fill_manual(values = c("#b2182b","#023858"))+
  labs(y="Normalized abundance(%)",x="20-TES-20");p
ggsave("../../result/TTS_result/metapicture.pdf",p , width = 9.24,height = 6.16)

