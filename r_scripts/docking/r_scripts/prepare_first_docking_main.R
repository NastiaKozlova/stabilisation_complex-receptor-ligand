part_name <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
#part_name<-part_name
part_scriprs<-strsplit(part_name,split ="/",fixed = T)[[1]]
part_scriprs<-paste(part_scriprs[1:(length(part_scriprs)-2)],collapse =  "/")
part_vina<-paste0(part_scriprs,"/programs/",collapse =  "")
part_scriprs<-paste0(part_scriprs,"/r_scripts/docking/r_scripts/")

part<-paste0(part_name,"docking_first/")
setwd(part)
v_ligand<-list.files(paste0("ligand/"))
a<-c()
for (i in 1:length(v_ligand)) {
  b<-strsplit(v_ligand[i],split = ".",fixed = T)[[1]][1]
  a<-c(a,b)
}
v_ligand<-a

df_ligand<-data.frame(matrix(ncol=2,nrow=length(v_ligand)))
colnames(df_ligand)<-c("ligand","c")
df_ligand$ligand<-v_ligand

df_active_center<-read.csv(paste0("active_center_surf.csv"),stringsAsFactors = F)
v_center<-unique(df_active_center$type)
df_center<-data.frame(matrix(ncol=2,nrow=length(v_center)))
colnames(df_center)<-c("center","c")
df_center$center<-v_center

df_ligand_center<-full_join(df_ligand,df_center,by="c",relationship = "many-to-many")
df_ligand_center$c<-NULL
write.csv(df_ligand_center,paste0("ligand_center_surf.csv"),row.names = F)

#docking
df_ligand_center<-read.csv(paste0("ligand_center_surf.csv"),stringsAsFactors = F)
df_ligand_center<-df_ligand_center%>%mutate(c="C")
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
df_all<-full_join(df_receptor,df_ligand_center,by="c",
                  relationship = "many-to-many")
df_all$c<-NULL
write.csv(df_all,paste0(part,"df_all_surf.csv"),row.names = F)
if (!dir.exists(paste0(part,"ligand/"))){dir.create(paste0(part,"ligand/"))}
if (!dir.exists(paste0(part,"receptor/"))){dir.create(paste0(part,"receptor/"))}
i<-1
for (i in 1:nrow(df_receptor)) {
  system(command = paste0(part_vina,"MGLTools-1.5.7/bin/pythonsh ",part_vina,"MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ",
                          part,"receptor_start/",df_receptor$receptor[i],".pdb -o ",part,"receptor/",df_receptor$receptor[i],".pdbqt ",
                          "-A None"))
}

system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_script.R ",part),ignore.stdout=T,wait = T)
