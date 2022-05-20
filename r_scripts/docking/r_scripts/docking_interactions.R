part_name <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-4

setwd(part_name)
setwd("din")


df_all<-read.csv(paste0("df_log_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(receptor_ligand=paste0(receptor,"_",ligand,"_",center))
if (!dir.exists(paste0("interaction/"))) { dir.create(paste0("interaction/"))}

i<-1
j<-1
for (j in 1:nrow(df_all)) {
  if(dir.exists(paste0("pdb_second/",df_all$receptor_ligand[j]))){
    if(!file.exists(paste0("interaction/",df_all$receptor_ligand[j],"/",df_all$new_number[j],".csv"))){
      a<-read.pdb(paste0(part_name,"receptor_start/",df_all$receptor[j],".pdb"))
      b<-read.pdb(paste0("pdb_second/",df_all$receptor_ligand[j],"/frame_",df_all$new_number[j],".pdb"))
      bs<-binding.site(a,b,cutoff=12)
      m<-bs$resnames
      a<-c()
      b<-c()
      y<-1
      for (y in 1:length(m)) {
        p<-strsplit(m[y],split = " ",fixed = T)[[1]][2]
        a<-c(a,p)
        p<-strsplit(m[y],split = " ",fixed = T)[[1]][1]
        b<-c(b,p)
      }
      a<-as.numeric(a)
      df_protein<-data.frame(matrix(ncol=2,nrow=length(a)))
      colnames(df_protein)<-c("resid","resno")
      df_protein$resid<-a
      df_protein$resno<-b
      if (!dir.exists(paste0("interaction/",df_all$receptor_ligand[j]))) { dir.create(paste0("interaction/",df_all$receptor_ligand[j]))}
      write.csv(df_protein,
                paste0("interaction/",df_all$receptor_ligand[j],"/",df_all$new_number[j],".csv"),
                row.names = F)
    }
  }
}