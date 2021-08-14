part_start <- commandArgs(trailingOnly=TRUE)
setwd(part_start)
library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)
v_rmsd<-4

v_first_bond<-list.files(paste0("din/interaction_fin/"))
if(length(v_first_bond)>0){
  df_first_bond_start<-read.csv(paste0("din/interaction_fin/",v_first_bond[1]),stringsAsFactors = F)
  for (i in 2:length(v_first_bond)) {
    df_first_bond_add<-read.csv(paste0("din/interaction_fin/",v_first_bond[i]),stringsAsFactors = F)
    df_first_bond_start<-rbind(df_first_bond_start,df_first_bond_add)
    df_first_bond_start<-df_first_bond_start%>%filter(persent_interactions>80)
  }
  df_first_bond_start<-df_first_bond_start%>%filter(persent_interactions>90)
#  df_first_bond_start<-df_first_bond_start%>%filter(persent_interactions==100)
  
  df_first_bond_start<-df_first_bond_start%>%mutate(receptor=NA)
  for (i in 1:nrow(df_first_bond_start)) {
    df_first_bond_start$receptor[i]<-strsplit(df_first_bond_start$system[i],split = "_",fixed = T)[[1]][1]  
  }
  df_first_bond_start<-df_first_bond_start%>%select(receptor,ligand, center, resid,resno, grops,grops_number,system)
  df_first_bond<-unique(df_first_bond_start)
  df_first_bond<-df_first_bond%>%filter(grops_number==1)  
  
  df_first_bond<-df_first_bond%>%mutate(test_complex_amino=paste(receptor,ligand, center, resid,resno,sep="_"))
  df_first_bond<-df_first_bond%>%group_by(test_complex_amino)%>%mutate(number_interactons_complex=n())
  df_first_bond<-ungroup(df_first_bond)

  df_first_bond_system<-df_first_bond%>%mutate(test_complex=paste(receptor,ligand, center, sep="_"))
  df_first_bond_system<-df_first_bond_system%>%select(receptor,ligand,center,test_complex,system,grops)
  df_first_bond_system<-unique(df_first_bond_system)
  
  df_first_bond_system<-df_first_bond_system%>%group_by(test_complex)%>%mutate(number_models_complex=n())
  df_first_bond_system<-ungroup(df_first_bond_system)
  df_first_bond_system<-df_first_bond_system%>%select(receptor, ligand,center,number_models_complex)
  df_first_bond_system<-unique(df_first_bond_system)
  
  df_first_bond<-ungroup(df_first_bond)
  df_first_bond_system<-ungroup(df_first_bond_system)
  df_first_bond_system$system<-NULL
  df_first_bond$system<-NULL

  df_first_bonda<-left_join(df_first_bond,df_first_bond_system,
                            by = c("receptor", "ligand", "center"))
  
  df_first_bonda<-df_first_bonda%>%mutate(occurence_interctions_complex=number_interactons_complex/number_models_complex*100)
  
  df_first_bonda<-df_first_bonda%>%select(receptor, ligand, center, resid, resno, occurence_interctions_complex)
  df_first_bonda<-unique(df_first_bonda)
  write.csv(df_first_bonda,"docking_aminoacid_interactions.csv",row.names = F)
}