setwd("/Users/zhangbing/work/zhaoslab/yulei/upload_code/script/polyU_analysis")
a=dir("../../data/polyU_data/")

for (each in a){
    
    command=sprintf(" Rscript  polyU.R  -f  \"../../data/polyU_data/%s/%s genome.fasta\"  -g \"../../data/polyU_data/%s/sequence.gff3\"  -o  '../../result/polyU_result/%s'  " , each,each,each,each)
    system(command  )
   cat(command)
    
}
