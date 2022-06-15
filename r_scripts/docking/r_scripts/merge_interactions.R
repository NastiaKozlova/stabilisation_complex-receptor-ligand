part_analysis <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-4

setwd(part_analysis)
setwd("din")

v_merged<-list.files("fin_merged_center/")
df_merged<-read.csv(paste0("fin_merged_center/",v_merged[1]),stringsAsFactors = F)
for (i in 2:length(v_merged)) {
  df_add<-read.csv(paste0("fin_merged_center/",v_merged[i]),stringsAsFactors = F)
  df_merged<-rbind(df_merged,df_add)
}

df_groups<-df_groups%>%mutate(name=paste0(ligand_center,"_",grop_number,"_",models.y))
df_merged<-df_merged%>%group_by(ligand_center)%>%mutate(number=n())


v_groups<-list.files("groups_fin/")
df_groups<-read.csv(paste0("groups_fin/",v_groups[1]),stringsAsFactors = F)
for (i in 2:length(v_groups)) {
  df_add<-read.csv(paste0("groups_fin/",v_groups[i]),stringsAsFactors = F)
  df_groups<-rbind(df_groups,df_add)
}
df_groups<-df_groups%>%mutate(name=paste0(ligand_center,"_",grop_number,"_",models.y))
df_all<-left_join(x = df_merged,y = df_groups,by=c("name.y"="name"))

#colnames(df_all)

#df_all<-unique(df_all)
df_log<-read.csv(paste0("log_fin.csv"),stringsAsFactors = F)

df_log<-df_log%>%select(models.x,models.y,grop_number,group,ligand_center,affinity,     
                        receptor,ligand,center,name)
#df_loga<-df_log[df_log$name%in%df_all$name.y,]

df_all<-df_all%>%mutate(receptor_ligand=paste0(receptor,"_",ligand,"_",center))

#if (dir.exists(paste0("interaction/"))) { system(command = paste0("rm -r ",part_analysis,"din/interaction/"))}
#if (dir.exists(paste0("interaction_TEMP/"))) {system(command = paste0("rm -r ",part_analysis,"din/interaction_TEMP/"))}
if (dir.exists(paste0("interaction_serf/"))) {system(command = paste0("rm -r ",part_analysis,"din/interaction_serf/"))}

if (!dir.exists(paste0("interaction_serf/"))) { dir.create(paste0("interaction_serf/"))}
i<-1
j<-3
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
              paste0("interaction_serf/",v_structure[j],".csv"),row.names = F)
  }
}