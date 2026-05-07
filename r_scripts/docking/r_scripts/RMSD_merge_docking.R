part_analysis <- commandArgs(trailingOnly=TRUE)
#group ligand structures
#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-3.5

print(Sys.time())

setwd(part_analysis)
df_all<-read.csv(paste0("df_all.csv"),stringsAsFactors = F)
part<-paste0(part_analysis,"din/")

if(dir.exists(paste0(part,"fin_merged"))) {system(command = paste0("rm -r ",part,"fin_merged"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"structure_merged"))) {system(command = paste0("rm -r ",part,"structure_merged"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"groups_merged"))) {system(command = paste0("rm -r ",part,"groups_merged"),ignore.stdout=T,wait = T)}
if (!dir.exists(paste0(part,"RMSD_merged"))) {dir.create(paste0(part,"RMSD_merged"))}
if (!dir.exists(paste0(part,"groups_merged"))) {dir.create(paste0(part,"groups_merged"))}
if (!dir.exists(paste0(part,"structure_merged"))) {dir.create(paste0(part,"structure_merged"))}
if (!dir.exists(paste0(part,"fin_merged"))) {dir.create(paste0(part,"fin_merged"))}

if(file.exists("din/log_fin.csv")){
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
  df_all<-df_all%>%mutate(test=!is.na(x))
  v_test<-  unique(df_all$test) 
  df_all$test<-NULL
  if(length(v_test)==1){
    df_merge<-left_join(df_all,df_all,by=c("receptor","ligand","searching_field"),relationship = "many-to-many")
    if(v_test){
      df_merge<-df_merge%>%filter(x.x<=x.y)
      df_merge<-df_merge%>%filter(y.x<=y.y)
      df_merge<-df_merge%>%filter(z.x<=z.y)
      
      df_merge<-df_merge%>%filter((x.x+searching_field)>=x.y)
      df_merge<-df_merge%>%filter((y.x+searching_field)>=y.y)
      df_merge<-df_merge%>%filter((z.x+searching_field)>=z.y)
      #df_structure_merge_start<-rbind(df_structure_merge_start,df_structure_merge)
    }
    df_merge<-df_merge%>%select(receptor,ligand,searching_field,center.x,center.y)
    setwd(part)
    df_log_fin<-read.csv("log_fin.csv",stringsAsFactors = F)
    df_log_fin<-df_log_fin%>%mutate(name=paste0(ligand_center,"_",grop_number,"_",models.x))
    
    
    df_RMSD<-df_log_fin%>%select(name,receptor,ligand,searching_field,center)
    df_RMSD<-unique(df_RMSD)
    df_calculation<-left_join(df_merge,df_RMSD,by=c("receptor","ligand","searching_field","center.x"="center"),
                              relationship = "many-to-many")
    df_calculation<-left_join(df_calculation,df_RMSD,by=c("receptor","ligand","searching_field","center.y"="center"),
                              relationship = "many-to-many")
    
    df_separate<-df_calculation%>%select(receptor,ligand)
    df_separate<-unique(df_separate)
    p<-1
    for (p in 1:nrow(df_separate)) {
      df_calculation_sep<-df_calculation%>%filter(receptor==df_separate$receptor[p])%>%
        filter(ligand==df_separate$ligand[p])
      df_calculation_sep<-df_calculation_sep%>%mutate(RMSD=NA)
      for (j in 1:nrow(df_calculation_sep)) {
        pdb_1<-read.pdb(paste0("str_fin/",df_calculation_sep$name.x[j]))
        pdb_2<-read.pdb(paste0("str_fin/",df_calculation_sep$name.y[j]))
        
        df_calculation_sep$RMSD[j]<-rmsd(pdb_1,pdb_2)
      }
      print(Sys.time())
      write.csv(df_calculation_sep,paste0("RMSD_merged/",df_separate$receptor[p],"_",df_separate$ligand[p],".csv"),row.names=F)
    }
    print(Sys.time())
  }
}