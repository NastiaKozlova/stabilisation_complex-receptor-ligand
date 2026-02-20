part_name <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
#part_name<-part_name
part_scriprs<-strsplit(part_name,split ="/",fixed = T)[[1]]
part_scriprs<-paste(part_scriprs[1:(length(part_scriprs)-1)],collapse =  "/")
part_vina<-paste0(part_scriprs,"/programs/",collapse =  "")
#part_scriprs<-paste0(part_scriprs,"/r_scripts/docking/r_scripts/")

v_part<-paste0(part_name,"/docking/docking_first/")
v_temp<-c("center","surf")
q<-1
for(q in 1:length(v_temp)){
  part<-paste0(v_part,v_temp[q],"/")
  setwd(part)
  df_ligand<-read.csv(paste0("ligand_field.csv"),stringsAsFactors = F)
  df_ligand<-df_ligand%>%mutate(c=NA)
  
  
  
  df_active_center<-read.csv(paste0("active_center.csv"),stringsAsFactors = F)
  df_center<-df_active_center%>%select(type)
  df_center<-unique(df_center)
  df_center<-df_center%>%mutate(c=NA)
  
  df_ligand_center<-full_join(df_ligand,df_center,by="c",relationship = "many-to-many")
  df_ligand_center$c<-NULL
  df_ligand_center<-unique(df_ligand_center)
  write.csv(df_ligand_center,paste0("ligand_center.csv"),row.names = F)
  
  #docking
  df_ligand_center<-read.csv(paste0("ligand_center.csv"),stringsAsFactors = F)
  df_ligand_center<-df_ligand_center%>%mutate(c="C")
  df_ligand_center<-df_ligand_center%>%mutate(center=type)
  v_receptor<-list.files(paste0("receptor_start/"))
  a<-c()
  for (i in 1:length(v_receptor)) {
    b<-strsplit(v_receptor[i],split = ".",fixed = T)[[1]][1]
    a<-c(a,b)
  }
  v_receptor<-a
  
  df_receptor<-data.frame(matrix(ncol=2,nrow=length(v_receptor)))
  colnames(df_receptor)<-c("receptor","c")
  df_receptor$receptor<-v_receptor
  df_receptor<-df_receptor%>%mutate(c="C")
  df_all<-full_join(df_receptor,df_ligand_center,by="c")
  df_all$c<-NULL
  write.csv(df_all,paste0("df_all.csv"),row.names = F)
  #if (!dir.exists(paste0("ligand/"))){dir.create(paste0("ligand/"))}
  #if (!dir.exists(paste0("receptor/"))){dir.create(paste0("receptor/"))}
  i<-1
  for (i in 1:nrow(df_receptor)) {
    system(command = paste0(part_vina,"MGLTools-1.5.7/bin/pythonsh ",part_vina,"MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ",
                            part,"receptor_start/",df_receptor$receptor[i],".pdb -o ",part,"receptor/",df_receptor$receptor[i],".pdbqt ",
                            "-A None"))
  }
  
  #system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_script.R ",part),ignore.stdout=T,wait = T)
}