part_start <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
setwd(part_start)
v_parta<-list.files("MD")
v_RMSD<-5
v_part<-paste0(part_start,"MD/",v_parta)
#part<-v_part[1]
num_model<-1
i<-1
j<-1

if (!dir.exists(paste0("MD_analysis/docking/"))){dir.create(paste0("MD_analysis/docking/"))}
if (!dir.exists(paste0("MD_analysis/docking/receptor"))){dir.create(paste0("MD_analysis/docking/receptor"))}
if (!dir.exists(paste0("MD_analysis/docking/receptor_start"))){dir.create(paste0("MD_analysis/docking/receptor_start"))}
if (!dir.exists(paste0("MD_analysis/docking/receptor_test"))){dir.create(paste0("MD_analysis/docking/receptor_test"))}
if (!dir.exists(paste0("MD_analysis/docking/df_RMSD"))){dir.create(paste0("MD_analysis/docking/df_RMSD"))}
max_num<-10
main_part<-c(8)
df_main<-data.frame(matrix(ncol=2,nrow = length(main_part)))
colnames(df_main)<-c("number","frames")
df_main$number<-main_part
for (j in 1:length(v_parta)) {
  part<-paste0(part_start,"MD/",v_parta[j])
  if (!dir.exists(paste0(part,"/din/docking/"))){dir.create(paste0(part,"/din/docking/"))}
  if (!dir.exists(paste0("MD_analysis/docking/receptor_test/",v_parta[j]))){dir.create(paste0("MD_analysis/docking/receptor_test/",v_parta[j]))}

  for (q in 1:nrow(df_main)) {
    df_main$frames[q]<-length(list.files(paste0(part,"/din/pdb_second/",df_main$number[q],"/")))-1
  }
  df_main<-df_main%>%filter(frames>0)
  df_main<-df_main%>%filter(number==max(number))
  df_ramachadran<-read.csv(paste0("MD_analysis/din/",v_parta[j],"/",df_main[1],"_time_Ramachadran.csv"),stringsAsFactors = F)
  df_ramachadran<-df_ramachadran%>%filter(ramachadran==min(ramachadran))
  for (i in 1:nrow(df_ramachadran)) {
    pdb<-read.pdb(paste0(part,"/din/pdb_second/",df_main$number[1],"/",df_ramachadran$frame_number[i],".pdb"))
    pdb.int<-atom.select(pdb,"water","ions",operator = "OR",inverse=T)
    pdb_fin<-trim.pdb(pdb,pdb.int)
    write.pdb(pdb_fin,paste0("MD_analysis/docking/receptor_test/",v_parta[j],"/",df_ramachadran$frame_number[i],".pdb"))
  }
  v_model<-list.files(paste0("MD_analysis/docking/receptor_test/",v_parta[j],"/"))
  df_RMSD<-data.frame(matrix(ncol = 2,nrow = length(v_model)))
  colnames(df_RMSD)<-c("model","RMSD")
  df_RMSD$model<-v_model
  df_RMSD<-full_join(df_RMSD,df_RMSD,by="RMSD")
  for (i in 1:nrow(df_RMSD)) {
    pdb_1<-read.pdb(paste0("MD_analysis/docking/receptor_test/",v_parta[j],"/",df_RMSD$model.x[i]))
    pdb_2<-read.pdb(paste0("MD_analysis/docking/receptor_test/",v_parta[j],"/",df_RMSD$model.y[i]))
    pdb_1.int<-atom.select(pdb_1,elety="CA")
    pdb_2.int<-atom.select(pdb_2,elety="CA")
#    df_RMSD$RMSD[i]<-rmsd(a = pdb_1,b = pdb_2,a.inds = pdb_1.int,b.inds = pdb_2.int,fit = T)
    df_RMSD$RMSD[i]<-rmsd(a = pdb_1,b = pdb_2,fit = T)
  }
  df_RMSD<-df_RMSD%>%filter(model.x!=model.y)
  df_RMSD<-df_RMSD%>%filter(RMSD<v_RMSD)
  df_RMSD<-df_RMSD%>%filter(RMSD<median(df_RMSD$RMSD))
  df_RMSD<-df_RMSD%>%group_by(model.x)%>%mutate(number_model=n())
  df_RMSD<-ungroup(df_RMSD)
  df_RMSD<-df_RMSD%>%filter(number_model>1)
  df_RMSD<-df_RMSD%>%select(model.x,number_model)
  df_RMSD<-unique(df_RMSD)
  df_RMSD<-df_RMSD%>%filter(number_model>(length(v_model)/2))
  df_RMSD<-df_RMSD%>%mutate(fin_model=paste0(v_parta[j],"_",model.x))
  for (i in 1:nrow(df_RMSD)) {
    df_RMSD$fin_model[i]<-strsplit( df_RMSD$fin_model[i],split = '.',fixed = T)[[1]][1]
  }
  write.csv(df_RMSD,paste0("MD_analysis/docking/df_RMSD/",v_parta[j],".csv"),row.names=F)
  for (i in 1:nrow(df_RMSD)) {
    pdb_1<-read.pdb(paste0("MD_analysis/docking/receptor_test/",v_parta[j],"/",df_RMSD$model.x[i]))
    write.pdb(pdb_1,paste0("MD_analysis/docking/receptor_start/",v_parta[j],"_",df_RMSD$model.x[i]))
  }
}
df_RMSD<-read.csv(paste0("MD_analysis/docking/df_RMSD/",v_parta[1],".csv"),stringsAsFactors = F)
if(length(v_parta)>1){
  for (j in 2:length(v_parta)) {
    df_RMSD_add<-read.csv(paste0("MD_analysis/docking/df_RMSD/",v_parta[j],".csv"),stringsAsFactors = F)
    df_RMSD<-rbind(df_RMSD,df_RMSD_add)
  }
}
write.csv(df_RMSD,paste0("MD_analysis/docking/df_RMSD.csv"),row.names=F)