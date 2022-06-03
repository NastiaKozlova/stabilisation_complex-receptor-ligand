part_analysis <- commandArgs(trailingOnly=TRUE)

#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-3.5

print(Sys.time())

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
df_all<-df_all%>%filter(!is.na(x))
df_all$name<-NULL

part<-paste0(part_analysis,"din/")
setwd(part)
if(dir.exists(paste0(part,"fin_merged"))) {system(command = paste0("rm -r ",part,"fin_merged"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"structure_merged"))) {system(command = paste0("rm -r ",part,"structure_merged"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"groups_merged"))) {system(command = paste0("rm -r ",part,"groups_merged"),ignore.stdout=T,wait = T)}
if (!dir.exists("RMSD_merged")) {dir.create("RMSD_merged")}
if (!dir.exists("groups_merged")) {dir.create("groups_merged")}
if (!dir.exists("structure_merged")) {dir.create("structure_merged")}
if (!dir.exists("fin_merged")) {dir.create("fin_merged")}


v_structure<-list.files(paste0(part_analysis,"din/str_fin"))
df_structure<-data.frame(matrix(nrow=length(v_structure),ncol = 5))
colnames(df_structure)<-c("name","receptor","ligand","center","RMSD")
df_structure$name<-v_structure
for (j in 1:nrow(df_structure)) {
  df_structure$receptor[j]<-strsplit(x = df_structure$name[j],split = "_")[[1]][1]
  df_structure$ligand[j]<-strsplit(x = df_structure$name[j],split = "_")[[1]][2]
}
v_center<-unique(df_all$center)
i<-1
for (i in 1:length(v_center)) {
  df_structure$center[grepl(x = df_structure$name,pattern = v_center[i])]<-v_center[i]
}
df_structure<-df_structure%>%filter(!is.na(center))

df_structure<-left_join(df_structure,df_all,by=c("receptor", "ligand", "center"))

df_analysis<-df_structure%>%select(receptor,ligand)
df_analysis<-unique(df_analysis)
df_analysis<-df_analysis%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
q<-1
i<-1
print(Sys.time())

for (q in 1:nrow(df_analysis)) {
  if(!file.exists(paste0("RMSD_merged/",df_analysis$receptor_ligand[q],".csv"))){
    df_structure_TEMP<-df_structure%>%filter(receptor==df_analysis$receptor[q])
    df_structure_TEMP<-df_structure_TEMP%>%filter(ligand==df_analysis$ligand[q])
    
    v_RMSD<-unique(df_structure_TEMP$center)
    df_structure_TEMPA<-df_structure_TEMP%>%filter(center==v_RMSD[1])
    df_structure_merge<-left_join(df_structure_TEMPA,df_structure_TEMP,by=c("receptor", "ligand","RMSD"))
    df_structure_merge_start<-df_structure_merge%>%filter(!is.na(RMSD))
    for (i in 1:length(v_RMSD)) {
      df_structure_TEMPA<-df_structure_TEMP%>%filter(center==v_RMSD[i])
      df_structure_merge<-left_join(df_structure_TEMPA,df_structure_TEMP,by=c("receptor", "ligand","RMSD"))
      df_structure_merge<-df_structure_merge%>%filter(x.x<=x.y)
      df_structure_merge<-df_structure_merge%>%filter(y.x<=y.y)
      df_structure_merge<-df_structure_merge%>%filter(z.x<=z.y)
      
      df_structure_merge<-df_structure_merge%>%filter((x.x+30)>=x.y)
      df_structure_merge<-df_structure_merge%>%filter((y.x+30)>=y.y)
      df_structure_merge<-df_structure_merge%>%filter((z.x+30)>=z.y)
      df_structure_merge_start<-rbind(df_structure_merge_start,df_structure_merge)
      
    }
    if(nrow(df_structure_merge_start)>0){
      for (j in 1:nrow(df_structure_merge_start)) {
        pdb_1<-read.pdb(paste0("str_fin/",df_structure_merge_start$name.x[j]))
        pdb_2<-read.pdb(paste0("str_fin/",df_structure_merge_start$name.y[j]))
        
        df_structure_merge_start$RMSD[j]<-rmsd(pdb_1,pdb_2)
      }
      df_structure_merge_start<-df_structure_merge_start%>%filter(RMSD<50)
    }
  
    write.csv(df_structure_merge_start,paste0("RMSD_merged/",df_analysis$receptor_ligand[q],".csv"),row.names=F)
    print(paste(df_analysis$receptor_ligand[q],Sys.time()))
  }
}
print(Sys.time())