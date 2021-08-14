part_start <- commandArgs(trailingOnly=TRUE)
#create new log file and save pdb_second
#print(Sys.time())
library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)
v_rmsd<-4

#part_start<-part_analysis
setwd(part_start)
df_all<-read.csv(paste0(part_start,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
num_model<-1
max_num<-5
if (!file.exists("din/")) {dir.create("din/")}
if (!dir.exists(paste0(part_start,"din/"))){dir.create(paste0(part_start,"din/"))}
if (!dir.exists(paste0(part_start,"din/log/"))){dir.create(paste0(part_start,"din/log/"))}
if (!dir.exists(paste0(part_start,"din/pdb_second/"))){dir.create(paste0(part_start,"din/pdb_second/"))}
a<-list.files(paste0("log/"))
df_topology<-data.frame(matrix(nrow=length(a),ncol =  3))
colnames(df_topology)<-c("name","run","name_log")

df_topology<-df_topology%>%mutate(exists="YES")


for (i in 1:length(a)) {
  b<-strsplit(a[i],split = ".",fixed = T)[[1]][1]
  
  df_topology$name_log[i]<-b
  b<-strsplit(b,split = "_")[[1]]
  df_topology$run[i]<-b[length(b)]
  df_topology$name[i]<-paste0(b[1:(length(b)-1)],collapse="_")
  

}

df_topology<-df_topology%>%filter(exists=="YES")
#print(nrow(df_topology))
df_topology<-left_join(df_topology,df_all,by="name")
df_topology<-df_topology%>%filter(!is.na(receptor))
write.csv(df_topology,"din/df_topology.csv",row.names = F)
print(Sys.time())
df_topology<-read.csv("din/df_topology.csv",stringsAsFactors = F)
df_ligand_center<-read.csv("ligand_center.csv",stringsAsFactors = F)
df_topology<-semi_join(df_topology,df_ligand_center,by=c("ligand", "center"))
print(Sys.time())
v_structure<-unique(df_topology$name)
for (j in 1:length(v_structure)) {
  if(!file.exists(paste0("din/log/",v_structure[j],".csv"))){
    df_topology_TEMP<-df_topology%>%filter(name==v_structure[j])
    df_log<-data.frame(matrix(ncol=5,nrow = 0))
    colnames(df_log)<-c("mode",   "affinity", "rmsd", "rmsd_from_BM","name_files")
    for (i in 1:nrow(df_topology_TEMP)) {
      if (file.size(paste0("log/",df_topology_TEMP$name_log[i],".log"))>1050){
        df_log_add<-read_table(paste0("log/",df_topology_TEMP$name_log[i],".log"),col_names = F,
                               skip = 25,n_max = 9,comment = "-",progress = F)
        if(ncol(df_log_add)>3){
          colnames(df_log_add)<-c("mode",   "affinity", "rmsd", "rmsd_from_BM")
          df_log_add<-df_log_add%>%mutate(name_files=df_topology_TEMP$name_log[i])
          df_log<-rbind(df_log,df_log_add)
        }
      }
    }
    df_log$mode<-as.numeric(df_log$mode)
    write.csv(df_log,paste0("din/log/",v_structure[j],".csv"),row.names = F)
  }
}
v_structure<-list.files("din/log/")
df_log<-read.csv(paste0("din/log/",v_structure[1]),stringsAsFactors = F)
for (j in 2:length(v_structure)) {

  df_log_add<-read.csv(paste0("din/log/",v_structure[j]),stringsAsFactors = F)
  if(nrow(df_log_add)>0){
    if(ncol(df_log_add)>3){
      df_log<-rbind(df_log,df_log_add)
    }
  }
}
print(Sys.time())
df_log$affinity<-as.numeric(df_log$affinity)
df_log<-df_log%>%filter(!is.na(affinity))
df_log<-left_join(df_log,df_topology,by=c("name_files"="name_log"))

df_log$mode<-as.numeric(df_log$mode)
print(Sys.time())
df_log<-df_log%>%group_by(name)%>%mutate(new_number=c(1:n()))
write.csv(df_log,"din/df_log_all.csv",row.names = F)
df_log<-read.csv("din/df_log_all.csv",stringsAsFactors = F)
if (!file.exists("din/pdb_second")) {dir.create("din/pdb_second")}
for (i in 1:nrow(df_log)){
  if (!file.exists(paste0("din/pdb_second/",df_log$name[i]))) {dir.create(paste0("din/pdb_second/",df_log$name[i]))}
  if (file.exists(paste0("analysis/",df_log$name_files[i],"_MODEL_",df_log$mode[i],".pdb"))){
    pdb<-read.pdb(paste0("analysis/",df_log$name_files[i],"_MODEL_",df_log$mode[i],".pdb"))
    pdb$atom$type<-"ATOM"
    write.pdb(pdb,paste0("din/pdb_second/",df_log$name[i],"/frame_",df_log$new_number[i],".pdb"))
  }
}
write.csv(df_log,"din/df_log_all.csv",row.names = F)
