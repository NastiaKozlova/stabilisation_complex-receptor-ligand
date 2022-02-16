part_start <- commandArgs(trailingOnly=TRUE)
#create new log file and save pdb_second

library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-4
#part_start<-part_analysis
setwd(part_start)
df_all<-read.csv(paste0(part_start,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
write.csv(df_all,"df_all.csv",row.names = F)
num_model<-1
max_num<-5
if (!file.exists("din/")) {dir.create("din/")}
if (!dir.exists(paste0(part_start,"din/"))){dir.create(paste0(part_start,"din/"))}
if (!dir.exists(paste0(part_start,"din/log/"))){dir.create(paste0(part_start,"din/log/"))}
if (!dir.exists(paste0(part_start,"din/pdb_second/"))){dir.create(paste0(part_start,"din/pdb_second/"))}
a<-list.files(paste0("out/"))
df_topology<-data.frame(matrix(nrow=length(a),ncol =  3))
colnames(df_topology)<-c("name","run","name_log")
i<-1
df_topology<-df_topology%>%mutate(exists="NO")
#i<-1
for (i in 1:length(a)) {
  b<-strsplit(a[i],split = ".",fixed = T)[[1]][1]
  df_topology$name_log[i]<-b
  b<-strsplit(b,split = "_")[[1]]
  df_topology$run[i]<-b[length(b)]
  df_topology$name[i]<-paste0(b[1:(length(b)-1)],collapse="_")
  if (file.exists(paste0("out/",df_topology$name_log[i],".pdbqt"))){
    df_topology$exists[i]<-"YES"
  }
}
df_topology<-df_topology%>%filter(!is.na(exists))
df_topology<-df_topology%>%filter(exists=="YES")
#print(nrow(df_topology))
df_topology<-left_join(df_topology,df_all,by="name")
df_topology<-df_topology%>%filter(!is.na(receptor))

df_log<-read.csv(paste0("din/log/",df_topology$name_log[1],".csv"),header = F)

colnames(df_log)<-c("mode",   "affinity", "rmsd", "rmsd_from_BM")
df_log<-df_log%>%mutate(name_files=df_topology$name_log[1])

for (i in 2:nrow(df_topology)) {
  df_log_add<-read.csv(paste0("din/log/",df_topology$name_log[i],".csv"),header = F)
  colnames(df_log_add)<-c("mode",   "affinity", "rmsd", "rmsd_from_BM")
  df_log_add<-df_log_add%>%mutate(name_files=df_topology$name_log[i])
  df_log<-rbind(df_log,df_log_add)
}
df_log$affinity<-as.numeric(df_log$affinity)
df_log<-df_log%>%filter(!is.na(affinity))
df_log<-left_join(df_log,df_topology,by=c("name_files"="name_log"))
#df_log<-df_log%>%mutate(ligand_center=paste0(ligand,"_",center,"_",system))

df_log$mode<-as.numeric(df_log$mode)

write.csv(df_log,"din/df_log_all.csv",row.names = F)

write.csv(df_topology,"din/df_topology.csv",row.names = F)
if (!file.exists("din/pdb_second")) {dir.create("din/pdb_second")}
v_analysis<-list.files("analysis")
df_analysis<-data.frame(matrix(ncol=3,nrow = length(v_analysis)))
colnames(df_analysis)<-c("files","name_files","mode")
df_analysis$files<-v_analysis
i<-1
for (i in 1:nrow(df_analysis)) {
  v<-strsplit(df_analysis$files[i],split = ".",fixed = T)[[1]]
  v<-strsplit(v,split = "_",fixed = T)[[1]]
  df_analysis$name_files[i]<-paste0(v[1:(length(v)-2)],collapse="_")
  df_analysis$mode[i]<-v[length(v)]
}
df_analysis<-df_analysis%>%mutate(mode=as.numeric(mode))
df_fin<-left_join(df_log,df_analysis,by=c("name_files","mode"))

df_fina<-df_fin%>%group_by(name)%>%mutate(new_number=1:n())
df_fina<-ungroup(df_fina)
i<-1
for (i in 1:nrow(df_fina)){
  if (!file.exists(paste0("din/pdb_second/",df_fina$name[i]))) {dir.create(paste0("din/pdb_second/",df_fina$name[i]))}
  pdb<-read.pdb(paste0("analysis/",df_fina$files[i]))
  write.pdb(pdb,paste0("din/pdb_second/",df_fina$name[i],"/frame_",df_fina$new_number[i],".pdb"))
}
write.csv(df_fina,"din/df_log_all.csv",row.names = F)