part_name <- commandArgs(trailingOnly=TRUE)
#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-3.5

setwd(part_name)
df_all<-read.csv(paste0(part_name,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))

#group_size<-round(nrow(df_all)/100,digits = 0)

part<-paste0(part_name,"din/")
setwd(part)
if(dir.exists(paste0(part,"fin_merged"))) {system(command = paste0("rm -r ",part,"fin_merged"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"structure_merged"))) {system(command = paste0("rm -r ",part,"structure_merged"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"groups_merged"))) {system(command = paste0("rm -r ",part,"groups_merged"),ignore.stdout=T,wait = T)}
if (!dir.exists("RMSD_merged")) {dir.create("RMSD_merged")}
if (!dir.exists("groups_merged")) {dir.create("groups_merged")}
if (!dir.exists("structure_merged")) {dir.create("structure_merged")}
if (!dir.exists("fin_merged")) {dir.create("fin_merged")}

v_protein_name<-unique(df_all$receptor)
v_structure_RMSD<-list.files(paste0("str_fin/"))
df_structure_RMSD<-data.frame(matrix(ncol=4,nrow = length(v_structure_RMSD)))
colnames(df_structure_RMSD)<-c("name","receptor","ligand","RMSD")
df_structure_RMSD$name<-v_structure_RMSD

for (j in 1:nrow(df_structure_RMSD)) {
  df_structure_RMSD$receptor[j]<-strsplit(x = df_structure_RMSD$name[j],split = "_")[[1]][1]
  df_structure_RMSD$ligand[j]<-strsplit(x = df_structure_RMSD$name[j],split = "_")[[1]][2]
}
df_analysis<-df_structure_RMSD%>%select(receptor,ligand)
df_analysis<-unique(df_analysis)
df_analysis<-df_analysis%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
q<-1
for (q in 1:nrow(df_analysis)) {
  df_structure_RMSD_TEMP<-df_structure_RMSD%>%filter(receptor==df_analysis$receptor[q])
  df_structure_RMSD_TEMP<-df_structure_RMSD_TEMP%>%filter(ligand==df_analysis$ligand[q])
  
  group_size<-round(nrow(df_structure_RMSD_TEMP)/100,digits = 0)
  if(!file.exists(paste0("RMSD_merged/",df_analysis$receptor_ligand[q],".csv"))){
    df_structure_RMSD_analysis<-left_join(df_structure_RMSD_TEMP,df_structure_RMSD_TEMP,by=c("receptor","ligand","RMSD"))
    df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(name.x!=name.y)
    
    for (j in 1:nrow(df_structure_RMSD_analysis)) {
      pdb_1<-read.pdb(paste0("str_fin/",df_structure_RMSD_analysis$name.x[j]))
      pdb_2<-read.pdb(paste0("str_fin/",df_structure_RMSD_analysis$name.y[j]))
      
      df_structure_RMSD_analysis$RMSD[j]<-rmsd(pdb_1,pdb_2)
    }
    write.csv(df_structure_RMSD_analysis,paste0("RMSD_merged/",df_analysis$receptor_ligand[q],".csv"),row.names=F)
  }
}
