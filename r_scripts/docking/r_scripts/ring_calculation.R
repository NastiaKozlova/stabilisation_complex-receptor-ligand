part_analysis = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
setwd(part_analysis)
part_programs<-strsplit(part_analysis,split = "/",fixed = T)[[1]]
part_programs<-paste0(part_programs[1:(length(part_programs)-3)],collapse = "/")
part_programs<-paste0(part_programs,"/")
#df_calculation<-read.csv(paste0(part_analysis,"start/df_calculation.csv"),stringsAsFactors =F)
i<-1
j<-1
#for (j in 1:nrow(df_calculation)) {
parta<-paste0(part_analysis)

#  part<-paste0(parta,df_calculation$complex_name[j])
#  setwd(part)
#  if(file.exists("patchdock/params.txt")){file.remove("patchdock/params.txt")}
models<-list.files("receptor_start")
df_RMSD<-data.frame(matrix(ncol = 2,nrow=length(models)))
colnames(df_RMSD)<-c("models","RMSD")
df_RMSD$models<-models

#  if (!dir.exists("structure/")){dir.create("structure/")}
if (!dir.exists("interactions/")){dir.create("interactions/")}
if (!dir.exists(paste0("ring2/"))){dir.create(paste0("ring2/"))}

df_RMSD<-df_RMSD%>%mutate(script=NA)
for (i in 1:nrow(df_RMSD)) {
  if (file.exists(paste0(part_analysis,"/receptor_start/",df_RMSD$models[i]))){
    if (!file.exists(paste0(part_analysis,"/ring2/",df_RMSD$models[i],".txt"))){
      df_RMSD$script[i]<-paste0("cd ", part_programs,"programs/dist/\n\n",
                                part_programs,"programs/dist/bin/Ring -i ",part_analysis,"receptor_start/",df_RMSD$models[i]," > ",part_analysis,"/ring2/",df_RMSD$models[i],".txt\n",
                                "rm ", part_analysis,"/structure_ring/",df_RMSD$models[i],"_fasta_A\n",
                                "rm ", part_analysis,"/structure_ring/",df_RMSD$models[i],"_modified\n")
    }
  }
}
df_RMSD<-df_RMSD%>%filter(!is.na(script))

write.csv(df_RMSD,paste0(part_analysis,"/script_","test",".txt"), row.names = F)
#df_RMSD<-df_RMSD%>%select(script)
df_RMSD[is.na(df_RMSD)]<-""
#  df_RMSD<-df_RMSD[1:round(nrow(df_RMSD)/2,digits = 0),]
#  df_RMSD<-df_RMSD[round(nrow(df_RMSD)/2,digits = 0):nrow(df_RMSD),]
#  df_RMSD<-df_RMSD[round(nrow(df_RMSD)/2,digits = 0):nrow(df_RMSD),]
#  df_RMSD<-df_RMSD[1:round(nrow(df_RMSD)/2,digits = 0),]
df_script<-df_RMSD%>%select(script)
write.table(df_script,paste0(part_analysis,"/script_","test",".txt"), quote=F,row.names = F,col.names = F,sep = "\n")

system(paste0("chmod +x ",part_analysis,"/script_","test",".txt"),ignore.stdout=T,wait = T)
system(paste0(part_analysis,"/script_","test",".txt"),ignore.stdout=T,wait = T)
