part_start <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
#part_start<-part_name
part_scriprs<-paste0(part_start,"r_scripts/")
part_start<-paste0(part_start,"docking_first/")
setwd(part_start)
if(!dir.exists(paste0(part_start,"analysis"))){dir.create(paste0(part_start,"analysis"))}
v_start<-list.files(paste0(part_start,"out"))
v_finish<-list.files(paste0(part_start,"analysis"))

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
for (i in 1:length(v_start)) {
  a<-strsplit(v_start[i],split = ".",fixed = T)[[1]][1]
  a<-paste0(part_start,"out/",a,".pdbqt ",part_start,"analysis/",a,"_MODEL_%d.pdb")
  system(command = paste0("python ", part_scriprs,"pdbqt_to_pdbs.py ",a),ignore.stdout=T,wait = T)
}
