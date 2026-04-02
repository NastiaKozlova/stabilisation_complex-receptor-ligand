part_analysis <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
#part_name<-part_name
#part_TEMP<-strsplit(part_name_zone,split = ",",fixed = T)[[1]]
#part_name<-part_TEMP[1]
#v_search_zone<-part_TEMP[2]
part_scriprs<-part_analysis
#part<-paste0(part_name,"docking_first/",v_search_zone,"/")
setwd(part_analysis)
if(!dir.exists(paste0("analysis"))){dir.create(paste0("analysis"))}
if(!dir.exists(paste0("din"))){dir.create(paste0("din"))}
if(!dir.exists(paste0("din/log"))){dir.create(paste0("din/log"))}
v_start<-list.files(paste0("out"))
v_finish<-list.files(paste0("analysis"))

b<-c()
i<-1
if (length(v_finish)>0){
  for (i in 1:length(v_finish)) {
    a<-strsplit(v_finish[i],split ="_MODEL_" )[[1]][1]
    b<-c(b,a)
    b<-unique(b)
  }
}
v_finish<-paste0(b,".pdbqt")
v_finish<-unique(v_finish)
v_start<-v_start[!v_start%in%v_finish]
if(length(v_start)>0){
  for (i in 1:length(v_start)) {
    a<-strsplit(v_start[i],split = ".",fixed = T)[[1]][1]
    a<-paste0(part_analysis,"out/",a,".pdbqt ",part_analysis,"analysis/",a,"_MODEL_%d.pdb")
    system(command = paste0("python3 ", part_scriprs,"pdbqt_to_pdbs.py ",a),ignore.stdout=T,wait = T)
  }
}
