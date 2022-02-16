part_analysis <- commandArgs(trailingOnly=TRUE)
#group ligand structures
part_start<-part_analysis
library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)
v_rmsd<-4
setwd(part_start)
part<-paste0(part_start,"din/")
setwd(part)

if (!dir.exists("groups")) {dir.create("groups")}

df_all<-read.csv(paste0(part_start,"df_all.csv"),stringsAsFactors = F)

if (!dir.exists("RMSD_analysis/")){dir.create("RMSD_analysis/")}
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
df_all<-df_all%>%mutate(number=0)
i<-1
print(Sys.time())
for (i in 1:nrow(df_all)) {
  models<-list.files(paste0("pdb_second/",df_all$name[i]))
  df_all$number[i]<-length(models)
}
print(Sys.time())
df_all<-df_all%>%filter(number>0)
v_str<-list.files(paste0("RMSD_analysis/"))
a<-c()
if(length(v_str)>0){
  for (i in 1:length(v_str)) {
    b<-strsplit(v_str[i],split = ".",fixed = T)[[1]][1]
    a<-c(a,b)
  }
  v_str<-a
  df_all<-df_all[!df_all$name%in%v_str,]
}

if(nrow(df_all)>0){
  for (i in 1:nrow(df_all)) {
    models<-list.files(paste0("pdb_second/",df_all$name[i]))
    if(length(models)>1){
      if(!file.exists(paste0("RMSD_analysis/",df_all$name[i],".csv"))){
        df_RMSD<-data.frame(matrix(ncol = 2,nrow=length(models)))
        colnames(df_RMSD)<-c("models","RMSD")
        df_RMSD$models<-models
        df_RMSD_all<-full_join(df_RMSD,df_RMSD,by="RMSD")
        for (j in 1:nrow(df_RMSD_all)) {
          pdb_1<-read.pdb(paste0("pdb_second/",df_all$name[i],"/",df_RMSD_all$models.x[j]))
          pdb_2<-read.pdb(paste0("pdb_second/",df_all$name[i],"/",df_RMSD_all$models.y[j]))
          df_RMSD_all$RMSD[j]<-rmsd(pdb_1,pdb_2)
        }
        write.csv(df_RMSD_all,paste0("RMSD_analysis/",df_all$name[i],".csv"),row.names = F)
      }
    }
  }
}

df_all<-read.csv(paste0(part_start,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
#sort to grops
i<-1
for (i in 1:nrow(df_all)) {

  if(file.exists(paste0("RMSD_analysis/",df_all$name[i],".csv"))){
    if (!dir.exists(paste0("groups/",df_all$name[i]))) {dir.create(paste0("groups/",df_all$name[i]))}
    df_RMSD_all<-read.csv(paste0("RMSD_analysis/",df_all$name[i],".csv"),stringsAsFactors = F)
    df_RMSD_all<-df_RMSD_all%>%filter(RMSD<4)
    df_RMSD_all<-df_RMSD_all%>%group_by(models.x)%>%mutate(number=n())
    df_RMSD_all<-ungroup(df_RMSD_all)                                             
    df_RMSD_all<-df_RMSD_all%>%filter(number>5)
    if (nrow(df_RMSD_all)>0){
      df_RMSD<-df_RMSD_all%>%select(models.x,number)
      df_RMSD<-unique(df_RMSD)
      df_RMSD<-df_RMSD%>%arrange(desc(number))
      df_RMSD_all<-df_RMSD_all
      for (j in 1:nrow(df_RMSD)) {
        if (!is.na(df_RMSD$models.x[j])) {
          df_RMSD_all_test<-df_RMSD_all%>%filter(models.x==df_RMSD$models.x[j])
          if (nrow(df_RMSD_all_test)>5) {
            df_RMSD_all_test<-df_RMSD_all_test%>%mutate(grop_number=j)
            write.csv(df_RMSD_all_test,paste0("groups/",df_all$name[i],"/grop_",j,".csv"),row.names = F) 
          }
        }
        
        df_RMSD_all$models.x[df_RMSD_all$models.x%in%c(df_RMSD_all_test$models.y,df_RMSD_all_test$models.x)]<-NA
        df_RMSD_all$models.x[df_RMSD_all$models.y%in%c(df_RMSD_all_test$models.y,df_RMSD_all_test$models.x)]<-NA
        df_RMSD_all<-df_RMSD_all%>%filter(!is.na(models.x))
        
        df_RMSD$models.x[df_RMSD$models.x%in%df_RMSD_all_test$models.y]<-NA
        df_RMSD$models.x[df_RMSD$models.x%in%df_RMSD_all_test$models.x]<-NA
      }
    }
  }
}


#combine all groups logs files
i<-3
if (!dir.exists("groups_fin")) {dir.create("groups_fin")}
for (i in 1:nrow(df_all)) {
  v_groups<-list.files(paste0("groups/",df_all$name[i]))
  if(length(v_groups)>0){
    df_groups<-read.csv(paste0("groups/",df_all$name[i],"/",v_groups[1]))
    df_groups<-df_groups%>%mutate(group=v_groups[1])
    df_groups<-df_groups%>%mutate(ligand_center=df_all$name[i])
    if (length(v_groups)>1) {
      for (j in 2:length(v_groups)) {
        df_groups_add<-read.csv(paste0("groups/",df_all$name[i],"/",v_groups[j]))
        df_groups_add<-df_groups_add%>%mutate(group=v_groups[j])
        df_groups_add<-df_groups_add%>%mutate(ligand_center=df_all$name[i])
        df_groups<-rbind(df_groups,df_groups_add)
      }
    }
    write.csv(df_groups,paste0("groups_fin/",df_all$name[i],".csv"),row.names = F)
  }
}

#copy pdb files to groups spb dir
if (!dir.exists("str")) {dir.create("str")}
if (!dir.exists("str_fin")) {dir.create("str_fin")}
i<-1
j<-1
k<-1
k<-2
#v_groups<-list.files(paste0("groups_fin/"))
for (j in 1:length(df_all$name)) {
  if(file.exists(paste0("groups_fin/",df_all$name[j],".csv"))){
    df_RMSD<-read.csv(paste0("groups_fin/",df_all$name[j],".csv"),stringsAsFactors = F)
    for (k in 1:nrow(df_RMSD)) {
      if (!dir.exists(paste0("str/",df_RMSD$ligand_center[k]))) { dir.create(paste0("str/",df_RMSD$ligand_center[k]))}
      if (!dir.exists(paste0("str/",df_RMSD$ligand_center[k],"/",df_RMSD$grop_number[k]))) {
        dir.create(paste0("str/",df_RMSD$ligand_center[k],"/",df_RMSD$grop_number[k]))}
      pdb<-read.pdb(paste0("pdb_second/",df_RMSD$ligand_center[k],"/",df_RMSD$models.y[k]))
      write.pdb(pdb,paste0("str/",df_RMSD$ligand_center[k],"/",df_RMSD$grop_number[k],"/",df_RMSD$models.y[k]))
    }
    df_RMSD<-df_RMSD%>%filter(models.y==models.x)
    for (q in 1:nrow(df_RMSD)){
        write.pdb(pdb,paste0("str_fin/",df_RMSD$ligand_center[q],"_",df_RMSD$grop_number[q],"_",df_RMSD$models.y[q]))
    }
  }
}


#energy bonding
i<-1
if (!dir.exists("log_fin")) {dir.create("log_fin")}
if (!dir.exists("plot")) {dir.create("plot")}
df_log<-read.csv("df_log_all.csv",stringsAsFactors = F)
df_log<-df_log%>%mutate(models.y=paste0("frame_",new_number,".pdb"))
for (i in 1:nrow(df_all)) {
  if(file.exists(paste0("groups_fin/",df_all$name[i],".csv"))){
    df_groups<-read.csv(paste0("groups_fin/",df_all$name[i],".csv"),stringsAsFactors = F)
    df_fin<-left_join(df_groups,df_log,by=c("models.y","ligand_center"="name"))
    write.csv(df_fin,paste0("log_fin/",df_all$name[i],".csv"),row.names = F)
  }
}
df_fin<-read.csv(paste0("log_fin/",df_all$name[1],".csv"),stringsAsFactors = F)
for (i in 2:length(df_all$name)) {
  if(file.exists(paste0("groups_fin/",df_all$name[i],".csv"))){
    df_fin_add<-read.csv(paste0("log_fin/",df_all$name[i],".csv"),stringsAsFactors = F)
    df_fin<-rbind(df_fin,df_fin_add)
  }
}
write.csv(df_fin,"log_fin.csv",row.names = F)
