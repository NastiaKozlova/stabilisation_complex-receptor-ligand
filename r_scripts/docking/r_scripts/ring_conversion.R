part_analysis = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
library(readr)
setwd(part_analysis)


# if(file.exists("receptor_start/params.txt")){file.remove("receptor_start/params.txt")}
models<-list.files("receptor_start")
df_RMSD<-data.frame(matrix(ncol = 2,nrow=length(models)))
colnames(df_RMSD)<-c("models","RMSD")
df_RMSD$models<-models
#  df_RMSD<-df_RMSD%>%mutate(chain=model_name[name])
if (!dir.exists("structure/")){dir.create("structure/")}
if (!dir.exists("interactions/")){dir.create("interactions/")}
if (!dir.exists(paste0("ring2/"))){dir.create(paste0("ring2/"))}
if (!dir.exists(paste0("ring/"))){dir.create(paste0("ring/"))}
if (!dir.exists(paste0("ring_mod/"))){dir.create(paste0("ring_mod/"))}
df_RMSD<-df_RMSD%>%mutate(name=paste0(models,".txt"))
v_rmsd<-list.files("ring2/")
df_RMSD<-df_RMSD[df_RMSD$name%in%v_rmsd,]
v_rmsd<-list.files("ring/")
df_RMSD<-df_RMSD[!df_RMSD$name%in%v_rmsd,]
if (nrow(df_RMSD)>0){
  for (i in 1:nrow(df_RMSD)) {
    if (!file.exists(paste0("ring/",df_RMSD$models[i],".txt"))){
      df_ring<-read_delim(paste0("ring2/",df_RMSD$models[i],".txt"), delim = "\t", skip = 11)
      df_ring<-df_ring[1:(which(df_ring$NodeId1%in%"NodeId")-1),]
      df_ring<-df_ring%>%select(NodeId1,Interaction,NodeId2,Distance,Angle,Energy)
      write.csv(df_ring,paste0("ring/",df_RMSD$models[i],".txt"),row.names = F)
    }
  }
}
i<-1


#  if(file.exists("receptor_start/params.txt")){file.remove("receptor_start/params.txt")}
models<-list.files("receptor_start")
df_RMSD<-data.frame(matrix(ncol = 2,nrow=length(models)))
colnames(df_RMSD)<-c("models","RMSD")
df_RMSD$models<-models
#  df_RMSD<-df_RMSD%>%mutate(chain=model_name[name])
if (!dir.exists("structure/")){dir.create("structure/")}
if (!dir.exists("interactions/")){dir.create("interactions/")}
if (!dir.exists(paste0("ring2/"))){dir.create(paste0("ring2/"))}
if (!dir.exists(paste0("ring/"))){dir.create(paste0("ring/"))}
if (!dir.exists(paste0("ring_mod/"))){dir.create(paste0("ring_mod/"))}
df_RMSD<-df_RMSD%>%mutate(name=paste0(models,".txt"))
df_RMSD<-data.frame(matrix(ncol = 2,nrow=length(models)))
colnames(df_RMSD)<-c("models","RMSD")
df_RMSD$models<-models
df_RMSD<-df_RMSD%>%mutate(name=paste0(models,".txt"))
v_rmsd<-list.files("ring2/")
for (i in 1:nrow(df_RMSD)) {
  if (!file.exists(paste0("ring_mod/",df_RMSD$models[i],".txt"))){
    if (file.exists(paste0("ring/",df_RMSD$models[i],".txt"))){
      #        pdb<-read.pdb(paste0("receptor.pdb"))
      #       df_pdb<-pdb$atom
      
      df_ring<-read.csv(paste0("ring/",df_RMSD$models[i],".txt"),stringsAsFactors = F)   
      df_ring<-df_ring%>%select(NodeId1,NodeId2,Interaction)
      df_ring<-unique(df_ring)
      #df_ring<-df_ring%>%mutate(chain_1="A")
      #df_ring<-df_ring%>%mutate(chain_2="A")
      for (j in 1:nrow(df_ring)) {
        v_separate_1<-strsplit(df_ring$NodeId1[j],split = ":",fixed = T)[[1]]
        v_separate_2<-strsplit(df_ring$NodeId2[j],split = ":",fixed = T)[[1]]
        df_ring$NodeId1[j]<-v_separate_1[2]
        df_ring$NodeId2[j]<-v_separate_2[2]
        df_ring$Interaction[j]<-strsplit(df_ring$Interaction[j],split = ":",fixed = T)[[1]][1]
      }
      df_ring<-df_ring%>%mutate(NodeId1=as.numeric(NodeId1))
      df_ring<-df_ring%>%mutate(NodeId2=as.numeric(NodeId2))
      write.csv(df_ring,paste0("ring_mod/",df_RMSD$models[i],".txt"),row.names = F)
    }
  }
}
