part_analysis <- commandArgs(trailingOnly=TRUE)
#group ligand structures
part_TEMP<-strsplit(part_analysis,split = ",")[[1]]
part_start<-part_TEMP[1]
v_rmsd<-as.numeric(part_TEMP[2])
library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)
library(ggplot2)

setwd(part_start)
part<-paste0(part_start,"din/")
if(file.exists(paste0(part_start,"din/df_log_all.csv"))){
  setwd(part)
  
  if(dir.exists(paste0(part,"groups_fin_merged"))) {system(command = paste0("rm -r ",part,"groups_fin_merged"),ignore.stdout=T,wait = T)}
  if(dir.exists(paste0(part,"groups_merged"))) {system(command = paste0("rm -r ",part,"groups_merged"),ignore.stdout=T,wait = T)}
  if(dir.exists(paste0(part,"str_merged"))) {system(command = paste0("rm -r ",part,"str_merged"),ignore.stdout=T,wait = T)}
  #if(dir.exists(paste0(part,"str_merge"))) {system(command = paste0("rm -r ",part,"str_merge"),ignore.stdout=T,wait = T)}
  if (!dir.exists("groups_fin_merged")) {dir.create("groups_fin_merged")}
  if (!dir.exists("groups_merged")) {dir.create("groups_merged")}
  if (!dir.exists("str_merged")) {dir.create("str_merged")}
  df_all<-read.csv(paste0(part_start,"df_all.csv"),stringsAsFactors = F)
  df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand))
  df_all<-df_all%>%select(name,receptor,ligand)
  df_all<-unique(df_all)
  #sort to grops
  for (i in 1:nrow(df_all)) {
    
    if(file.exists(paste0("RMSD_merged/",df_all$name[i],".csv"))){
      print(df_all$name[i])
      if (!dir.exists(paste0("groups_merged/",df_all$name[i]))) {dir.create(paste0("groups_merged/",df_all$name[i]))}
      df_RMSD_all<-read.csv(paste0("RMSD_merged/",df_all$name[i],".csv"),stringsAsFactors = F)
      df_RMSD_all<-df_RMSD_all%>%filter(RMSD<v_rmsd)
      df_RMSD_all<-df_RMSD_all%>%group_by(name.x)%>%mutate(number=n())
      df_RMSD_all<-ungroup(df_RMSD_all)   
      df_RMSD_all<-df_RMSD_all%>%filter(number>1)
      #df_RMSD_all<-df_RMSD_all%>%filter(number>median(df_RMSD_all$number))
      if (nrow(df_RMSD_all)>0){
        df_RMSD<-df_RMSD_all%>%select(name.x,number)
        df_RMSD<-unique(df_RMSD)
        df_RMSD<-df_RMSD%>%arrange(desc(number))
        df_RMSD_all<-df_RMSD_all
        for (j in 1:nrow(df_RMSD)) {
          if (!is.na(df_RMSD$name.x[j])) {
            df_RMSD_all_test<-df_RMSD_all%>%filter(name.x==df_RMSD$name.x[j])
            if (nrow(df_RMSD_all_test)>1) {
              df_RMSD_all_test<-df_RMSD_all_test%>%mutate(grop_number=j)
              write.csv(df_RMSD_all_test,paste0("groups_merged/",df_all$name[i],"/grop_",j,".csv"),row.names = F) 
            }
          }
          
          df_RMSD_all$name.x[df_RMSD_all$name.x%in%c(df_RMSD_all_test$name.y,df_RMSD_all_test$name.x)]<-NA
          df_RMSD_all$name.x[df_RMSD_all$name.y%in%c(df_RMSD_all_test$name.y,df_RMSD_all_test$name.x)]<-NA
          df_RMSD_all<-df_RMSD_all%>%filter(!is.na(name.x))
          
          df_RMSD$name.x[df_RMSD$name.x%in%df_RMSD_all_test$name.y]<-NA
          df_RMSD$name.x[df_RMSD$name.x%in%df_RMSD_all_test$name.x]<-NA
        }
      }
    }
  }
  
  
  #combine all groups logs files
  i<-1
  if (!dir.exists("groups_fin_merged")) {dir.create("groups_fin_merged")}
  for (i in 1:nrow(df_all)) {
    v_groups<-list.files(paste0("groups_merged/",df_all$name[i]))
    if(length(v_groups)>0){
      df_groups<-read.csv(paste0("groups_merged/",df_all$name[i],"/",v_groups[1]))
      df_groups<-df_groups%>%mutate(group=v_groups[1])
      df_groups<-df_groups%>%mutate(ligand_center=df_all$name[i])
      if (length(v_groups)>1) {
        for (j in 2:length(v_groups)) {
          df_groups_add<-read.csv(paste0("groups_merged/",df_all$name[i],"/",v_groups[j]))
          df_groups_add<-df_groups_add%>%mutate(group=v_groups[j])
          df_groups_add<-df_groups_add%>%mutate(ligand_center=df_all$name[i])
          df_groups<-rbind(df_groups,df_groups_add)
        }
      }
      write.csv(df_groups,paste0("groups_fin_merged/",df_all$name[i],".csv"),row.names = F)
    }
  }
  
  #copy pdb files to groups spb dir
  #if (!dir.exists("str")) {dir.create("str")}
  if (!dir.exists("str_merged")) {dir.create("str_merged")}
  i<-1
  j<-1
  k<-1
  k<-2
  #v_groups<-list.files(paste0("groups_fin/"))
  for (j in 1:nrow(df_all)) {
    if(file.exists(paste0("groups_fin_merged/",df_all$name[j],".csv"))){
      df_RMSD<-read.csv(paste0("groups_fin_merged/",df_all$name[j],".csv"),stringsAsFactors = F)
      df_RMSD<-df_RMSD%>%filter(name.y==name.x)
#      for (q in 1:nrow(df_RMSD)){
#        pdb<-read.pdb(paste0("str_fin/",df_RMSD$name.y[q]))
        
#        write.pdb(pdb,paste0("str_merged/",df_RMSD$name.y[q]))
#      }
    }
  }
  
  
  #energy bonding
  i<-1
  if (!dir.exists("log_merged")) {dir.create("log_merged")}
  if (!dir.exists("plot")) {dir.create("plot")}
  df_log<-read.csv("log_fin.csv",stringsAsFactors = F)
  df_log<-df_log%>%mutate(name=paste0(ligand_center,"_",grop_number,"_",models.x))
  for (i in 1:nrow(df_all)) {
    if(file.exists(paste0("groups_fin_merged/",df_all$name[i],".csv"))){
      df_groups<-read.csv(paste0("groups_fin_merged/",df_all$name[i],".csv"),stringsAsFactors = F)
      df_fin<-left_join(df_groups,df_log,by=c("receptor",        "ligand" ,         "searching_field","name.y"="name"))
      write.csv(df_fin,paste0("log_merged/",df_all$name[i],".csv"),row.names = F)
    }
  }
  df_fin<-read.csv(paste0("log_merged/",df_all$name[1],".csv"),stringsAsFactors = F)
  if(nrow(df_all)>1){
    for (i in 2:nrow(df_all)) {
      if(file.exists(paste0("log_merged/",df_all$name[i],".csv"))){
        df_fin_add<-read.csv(paste0("log_merged/",df_all$name[i],".csv"),stringsAsFactors = F)
        df_fin<-rbind(df_fin,df_fin_add)
      }
    }
  }
  write.csv(df_fin,"log_merge.csv",row.names = F)
#  p<-ggplot(data = df_fin)+
#    geom_density(aes(x=affinity))+
#    facet_grid(ligand~center)+
#    theme_bw()
}