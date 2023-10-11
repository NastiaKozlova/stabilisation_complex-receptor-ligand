part_analysis <- commandArgs(trailingOnly=TRUE)
#group ligand structures
part_TEMP<-strsplit(part_analysis,split = ",")[[1]]
part_start<-part_TEMP[1]
v_rmsd<-10
v_str<-0
library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)

setwd(part_start)
part<-paste0(part_start,"din/")
setwd(part)
if(dir.exists(paste0(part,"groups"))) {system(command = paste0("rm -r ",part,"groups"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"groups_fin"))) {system(command = paste0("rm -r ",part,"groups_fin"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"str"))) {system(command = paste0("rm -r ",part,"str"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"str_fin_center"))) {system(command = paste0("rm -r ",part,"str_fin_center"),ignore.stdout=T,wait = T)}
if (!dir.exists("groups")) {dir.create("groups")}
if (!dir.exists("groups_fin")) {dir.create("groups_fin")}
if (!dir.exists("str_fin_center")) {dir.create("str_fin_center")}
df_all<-read.csv(paste0(part_start,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
df_all<-df_all%>%mutate(stuture_number=NA)
#sort to grops
for (i in 1:nrow(df_all)) {
  
  if(file.exists(paste0("RMSD_analysis/",df_all$name[i],".csv"))){
    print(df_all$name[i])
    if (!dir.exists(paste0("groups/",df_all$name[i]))) {dir.create(paste0("groups/",df_all$name[i]))}
    df_RMSD_all<-read.csv(paste0("RMSD_analysis/",df_all$name[i],".csv"),stringsAsFactors = F)
    df_RMSD_all<-df_RMSD_all%>%filter(RMSD<v_rmsd)
    df_all$stuture_number[i]<-nrow(df_RMSD_all)
    df_RMSD_all<-df_RMSD_all%>%group_by(models.x)%>%mutate(number=n())
    df_RMSD_all<-ungroup(df_RMSD_all)                                             
    df_RMSD_all<-df_RMSD_all%>%filter(number>v_str)
    if (nrow(df_RMSD_all)>v_str){
      df_RMSD<-df_RMSD_all%>%select(models.x,number)
      df_RMSD<-unique(df_RMSD)
      df_RMSD<-df_RMSD%>%arrange(desc(number))
      df_RMSD_all<-df_RMSD_all
      for (j in 1:nrow(df_RMSD)) {
        if (!is.na(df_RMSD$models.x[j])) {
          df_RMSD_all_test<-df_RMSD_all%>%filter(models.x==df_RMSD$models.x[j])
          if (nrow(df_RMSD_all_test)>0) {
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
#ggplot(data=df_all)+
#geom_freqpoly(aes(x=stuture_number))

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
#if (!dir.exists("str")) {dir.create("str")}
if (!dir.exists("str_fin")) {dir.create("str_fin")}
i<-1
j<-1
k<-1
k<-2
#v_groups<-list.files(paste0("groups_fin/"))
for (j in 1:length(df_all$name)) {
  if(file.exists(paste0("groups_fin/",df_all$name[j],".csv"))){
    df_RMSD<-read.csv(paste0("groups_fin/",df_all$name[j],".csv"),stringsAsFactors = F)
    df_RMSD<-df_RMSD%>%filter(models.y==models.x)
    for (q in 1:nrow(df_RMSD)){
      pdb<-read.pdb(paste0("pdb_second/",df_RMSD$ligand_center[q],"/",df_RMSD$models.y[q]))
      
      write.pdb(pdb,paste0("str_fin_center/",df_RMSD$ligand_center[q],"_",df_RMSD$grop_number[q],"_",df_RMSD$models.y[q]))
    }
  }
}


#energy bonding
i<-1
if (!dir.exists("log_fin_center")) {dir.create("log_fin_center")}
if (!dir.exists("plot")) {dir.create("plot")}
df_log<-read.csv("df_log_all_center.csv",stringsAsFactors = F)
df_log<-df_log%>%mutate(models.y=paste0("frame_",new_number,".pdb"))
for (i in 1:nrow(df_all)) {
  if(file.exists(paste0("groups_fin/",df_all$name[i],".csv"))){
    df_groups<-read.csv(paste0("groups_fin/",df_all$name[i],".csv"),stringsAsFactors = F)
    df_fin<-left_join(df_groups,df_log,by=c("models.y","ligand_center"="name"))
    write.csv(df_fin,paste0("log_fin_center/",df_all$name[i],".csv"),row.names = F)
  }
}
a<-list.files("log_fin_center")
df_all<-df_all[paste0(df_all$name,".csv")%in%a,]

df_fin<-read.csv(paste0("log_fin_center/",df_all$name[1],".csv"),stringsAsFactors = F)
for (i in 2:length(df_all$name)) {
  if(file.exists(paste0("groups_fin/",df_all$name[i],".csv"))){
    df_fin_add<-read.csv(paste0("log_fin_center/",df_all$name[i],".csv"),stringsAsFactors = F)
    df_fin<-rbind(df_fin,df_fin_add)
  }
}
write.csv(df_fin,"log_fin_center.csv",row.names = F)
