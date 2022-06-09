part_analysis <- commandArgs(trailingOnly=TRUE)
#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)
#v_rmsd<-3.5

setwd(part_analysis)
df_all<-read.csv(paste0(part_analysis,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
df_all<-df_all%>%mutate(x=NA)
df_all<-df_all%>%mutate(y=NA)
df_all<-df_all%>%mutate(z=NA)
i<-1
for (i in 1:nrow(df_all)) {
  a<-strsplit(df_all$center[i],split = "_")[[1]]
  df_all$x[i]<-as.numeric(a[3])
  df_all$y[i]<-as.numeric(a[5])
  df_all$z[i]<-as.numeric(a[7])
}
df_all<-df_all%>%filter(is.na(x))
#df_all$name<-NULL
df_all<-df_all%>%select(name,receptor,ligand,center)
#group_size<-round(nrow(df_all)/100,digits = 0)

part<-paste0(part_analysis,"din/")
setwd(part)
if(dir.exists(paste0(part,"fin_merged_center"))) {system(command = paste0("rm -r ",part,"fin_merged_center"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"structure_merged_center"))) {system(command = paste0("rm -r ",part,"structure_merged_center"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"groups_merged_center"))) {system(command = paste0("rm -r ",part,"groups_merged_center"),ignore.stdout=T,wait = T)}
if (!dir.exists("RMSD_merged_center")) {dir.create("RMSD_merged_center")}
if (!dir.exists("groups_merged_center")) {dir.create("groups_merged_center")}
if (!dir.exists("structure_merged_center")) {dir.create("structure_merged_center")}
if (!dir.exists("fin_merged_center")) {dir.create("fin_merged_center")}

v_protein_name<-unique(df_all$receptor)
v_structure_RMSD<-list.files(paste0("str_fin/"))
df_structure_RMSD<-data.frame(matrix(ncol=5,nrow = length(v_structure_RMSD)))
colnames(df_structure_RMSD)<-c("name","receptor","ligand","center","RMSD")
df_structure_RMSD$name<-v_structure_RMSD
j<-1
a<-c()
for (j in 1:nrow(df_all)) {
  b<-df_structure_RMSD$name[grepl(x = df_structure_RMSD$name,pattern = df_all$name[j])]
  a<-c(a,b)
}
df_structure_RMSD<-df_structure_RMSD[df_structure_RMSD$name%in%a,]
for (j in 1:nrow(df_structure_RMSD)) {
  a<-strsplit(x = df_structure_RMSD$name[j],split = "_frame")[[1]][1]
  a<-strsplit(x = a,split = "_")[[1]]
  df_structure_RMSD$receptor[j]<-a[1]
  df_structure_RMSD$ligand[j]<-a[2]
  df_structure_RMSD$center[j]<-paste0(a[3:(length(a)-1)],collapse = "_")
}
df_structure_RMSD<-df_structure_RMSD[df_structure_RMSD$center%in%c(df_all$center),]

df_analysis<-df_structure_RMSD%>%select(receptor,ligand)
df_analysis<-unique(df_analysis)
df_analysis<-df_analysis%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
q<-1
for (q in 1:nrow(df_analysis)) {
  df_structure_RMSD_TEMP<-df_structure_RMSD%>%filter(receptor==df_analysis$receptor[q])
  df_structure_RMSD_TEMP<-df_structure_RMSD_TEMP%>%filter(ligand==df_analysis$ligand[q])
  if(!file.exists(paste0("RMSD_merged_center/",df_analysis$receptor_ligand[q],".csv"))){
    df_structure_RMSD_analysis<-left_join(df_structure_RMSD_TEMP,df_structure_RMSD_TEMP,by=c("receptor","ligand","RMSD"))
    for (j in 1:nrow(df_structure_RMSD_analysis)) {
      pdb_1<-read.pdb(paste0("str_fin/",df_structure_RMSD_analysis$name.x[j]))
      pdb_2<-read.pdb(paste0("str_fin/",df_structure_RMSD_analysis$name.y[j]))
      
      df_structure_RMSD_analysis$RMSD[j]<-rmsd(pdb_1,pdb_2)
    }
    write.csv(df_structure_RMSD_analysis,paste0("RMSD_merged_center/",df_analysis$receptor_ligand[q],".csv"),row.names=F)
  }
}
