part_start <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)
v_rmsd<-4

setwd(part_start)

main_part<-list.files("interaction")
if (!dir.exists("din/interaction_fin/")){dir.create("din/interaction_fin/")}
i<-1
j<-643
df_topology<-read.csv("din/df_topology.csv",stringsAsFactors = F)
v_sort_list<-list.files(paste0("din/interaction/"))
df_topology<-df_topology[df_topology$name%in%v_sort_list,]
df_topology<-df_topology%>%select(name,receptor,ligand,center)
df_topology<-unique(df_topology)
if(nrow(df_topology)>0){
  for (j in 1:nrow(df_topology)) {
    pdb<-read.pdb(paste0("receptor_start/",df_topology$receptor[j],".pdb"))
    v_groups<-list.files(paste0("din/interaction/",df_topology$name[j]))
    if(length(v_groups)>0){
      for (i in 1:length(v_groups)) {
        v_frame<-list.files(paste0("din/interaction/",df_topology$name[j],"/",v_groups[i]))
        df_pdb<-pdb$atom
        df_pdb<-df_pdb%>%filter(elety=="CA")
        df_pdb<-df_pdb%>%mutate(number_interactions=0)
        df_pdb<-df_pdb%>%mutate(system=df_topology$receptor[j])
        df_pdb<-df_pdb%>%mutate(center=df_topology$center[j])
        df_pdb<-df_pdb%>%mutate(ligand=df_topology$ligand[j])
        df_pdb<-df_pdb%>%mutate(grops=v_groups[i])
        df_pdb<-df_pdb%>%mutate(grops_number=i)
        for (q in 1:length(v_frame)) {
          df_interaction<-read.csv(paste0("din/interaction/",df_topology$name[j],"/",v_groups[i],"/",v_frame[q]),stringsAsFactors = F)
          colnames(df_interaction)<-c(colnames(df_interaction)[2],colnames(df_interaction)[1])
          df_pdb$number_interactions[df_pdb$resno%in%df_interaction$resno]<-df_pdb$number_interactions[df_pdb$resno%in%df_interaction$resno]+1
        }
        df_pdb<-df_pdb%>%mutate(persent_interactions=number_interactions/length(v_frame)*100)
        write.csv(df_pdb,paste0("din/interaction_fin/",df_topology$name[j],"_",v_groups[i],".csv"),row.names = F)
      }
    } 
  }
  v_groups<-list.files(paste0("din/interaction_fin/"))
  df_pdb<-read.csv(paste0("din/interaction_fin/",v_groups[1]),stringsAsFactors =  F)
  for (i in 2:length(v_groups)) {
    df_pdb_add<-read.csv(paste0("din/interaction_fin/",v_groups[i]),stringsAsFactors =  F)
    df_pdb<-rbind(df_pdb,df_pdb_add)
  }
  
  df_pdb<-df_pdb%>%mutate(aminoacids=paste(resno,resid))
  df_pdb<-df_pdb%>%select(resno, x, y, z, number_interactions, system, center, ligand, 
                          grops, grops_number,persent_interactions, aminoacids )
  write.csv(df_pdb,"interaction_fin.csv",row.names = F)
}