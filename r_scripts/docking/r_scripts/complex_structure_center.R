part_analysis <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(dplyr)
library(ggplot2)
paste0(part_analysis)

df_merge<-read.csv(paste0(part_analysis,"din/","df_merge_structure_log_center.csv"),stringsAsFactors = F)
#df_merge<-df_merge%>%mutate(name=paste0(ligand,"_",center,"_",models.x))
#df_merge<-semi_join(df_merge,df_all)
df_merge<-df_merge%>%select(name.x,receptor,ligand,size_of_group)
df_merge<-unique(df_merge)
i<-1
if(dir.exists(paste0(part_analysis,"din/complex_structure_center"))) {system(command = paste0("rm -r ",part_analysis,"din/complex_structure_center"),ignore.stdout=T,wait = T)}

if(!dir.exists(paste0(part_analysis,"din/complex_structure_center"))){dir.create(paste0(part_analysis,"din/complex_structure_center"))}
for (i in 1:nrow(df_merge)) {
  receptor_name<-paste0(part_analysis,"receptor_start/",df_merge$receptor[i],".pdb")
  ligand_name<-paste0(part_analysis,"din/str_fin/",df_merge$name.x[i])
  pdb_receptor<-read.pdb(receptor_name)
  pdb_ligand<-read.pdb(ligand_name)
  pdb_complex<-cat.pdb(pdb_receptor, pdb_ligand, rechain=TRUE)
  write.pdb(pdb_complex,paste0(part_analysis,"din/complex_structure_center/",df_merge$name.x[i]))
}