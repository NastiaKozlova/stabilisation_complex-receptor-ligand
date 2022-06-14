part_analysis <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-4

setwd(part_analysis)
setwd("din")


df_all<-read.csv(paste0("df_merge_structure_log_center.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(receptor_ligand=paste0(receptor,"_",ligand,"_",center))

if (dir.exists(paste0("interaction_center/"))) {system(command = paste0("rm -r ",part_analysis,"din/interaction_center/"))}

if (!dir.exists(paste0("interaction_center/"))) { dir.create(paste0("interaction_center/"))}
i<-1
j<-1
p<-1
v_structure<-unique(df_all$name.x)
for (j in 1:length(v_structure)) {
  df_complex<-df_all%>%filter(name.x==v_structure[j])
  pdb<-read.pdb(paste0(part_analysis,"receptor_start/",df_all$receptor[j],".pdb"))
  
  df_pdb<-pdb$atom
  df_pdb<-df_pdb%>%filter(elety=="CA")
  df_pdb<-df_pdb%>%mutate(number_interactions=0)
  df_pdb<-df_pdb%>%mutate(tested_structure=0)
  df_pdb<-df_pdb%>%mutate(total_structure=nrow(df_complex))
  test<-nrow(df_pdb)
  for (p in 1:nrow(df_complex)) {
    if(file.exists(paste0("interaction/",df_complex$receptor_ligand[p],"/",df_complex$new_number[p],".csv"))){
      df_protein<-read.csv(paste0("interaction/",df_complex$receptor_ligand[p],"/",df_complex$new_number[p],".csv"),
                           stringsAsFactors = F) 
      df_pdb$number_interactions[df_pdb$resno%in%df_protein$resid]<-df_pdb$number_interactions[df_pdb$resno%in%df_protein$resid]+1
      df_pdb$tested_structure<-df_pdb$tested_structure+1
      
    }
  }
  df_pdb<-df_pdb%>%filter(tested_structure==total_structure)
  if(nrow(df_pdb)==test){
    df_pdb<-df_pdb%>%select(resno,resid,x,y,z,number_interactions,tested_structure,total_structure)
    df_pdb<-df_pdb%>%mutate(persent_interactions=number_interactions/total_structure*100)
    write.csv(df_pdb,
              paste0("interaction_center/",v_structure[j],".csv"),row.names = F)
  }
}