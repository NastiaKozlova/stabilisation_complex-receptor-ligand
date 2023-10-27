part_prepare = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
part_TEMP<-strsplit(part_prepare,split = ",",fixed = T)[[1]]
part_name<-part_TEMP[1]
v_list_protein<-part_TEMP[2]
part_prepare<-paste0(part_name,v_list_protein,"/")
setwd(part_prepare)
v_ligands<-list.files(paste0("MD_globular_protein_ligand/start/ligands_center/"))
v_receptors<-list.files(paste0("MD_globular_protein_ligand/start/receptor"))
df_resname_change<-read.csv(paste0(part_name,"start/receptor_ligand_MD/resname_mutate.csv"),stringsAsFactors = F)
df_active_center<-read.csv(paste0("docking/docking_first/din/df_merge_structure_log_center.csv"),stringsAsFactors =  F)
df_active_center<-df_active_center%>%group_by(name.x)%>%mutate(group_size=n())
df_active_center<-df_active_center%>%select(name.x,receptor,ligand,center.x,
                                            group_size)
df_active_center<-ungroup(df_active_center)
df_active_center<-unique(df_active_center)
df_active_center<-df_active_center%>%mutate(name=paste0(receptor,"_",ligand,"_",center.x))
df_active_center<-df_active_center%>%filter(ligand!="chloramphenicol")
df_active_center<-df_active_center%>%filter(ligand!="norfloxacin")
#df_active_center<-left_join(df_active_center,df_resname_change,by=c("ligand"="ligand_name"))

#df_active_center<-df_active_center%>%filter(!is.na(charmm_name))
#df_active_center<-df_active_center%>%filter("models.x"      "RMSD"          "models.y"      "number"        "grop_number"   "group"         "ligand_center" "receptor"     
#                                            [9] "ligand"        "center"        "size_of_group" "name")
#<-data.frame(matrix(ncol = 4,nrow = length(v_ligands)))
#colnames(df_ligands)<-c("file_name","structure_name","complex_name","resname_change")
#df_ligands$file_name<-v_ligands
#i<-1
#for (i in 1:nrow(df_ligands)) {
#  a<-strsplit(df_ligands$file_name[i],split = "_")[[1]]
#  df_ligands$complex_name[i]<-paste0(a[1:(length(a)-3)],collapse =  "_")
#  df_ligands$structure_name[i]<-paste0(a[(length(a)-2):(length(a))],collapse =  "_")
#}
#df_ligands<-left_join(df_ligands,df_active_center,by = c("complex_name"="name"))
#i<-1
#for (i in 1:nrow(df_ligands)) {
#  pdb<-read.pdb(paste0("MD_globular_protein_ligand/start/ligands/",df_ligands$file_name[i]))
#  if(length(unique(pdb$atom$resid))==1){
#    df_ligands$resname_change<-unique(pdb$atom$resid)
#  }
#}
#df_ligands<-df_ligands%>%filter(!is.na(resname_change))
#df_ligands<-df_ligands%>%filter(!is.na(center))
#df_ligands<-df_ligands%>%group_by(ligand)%>%mutate(number=1:n())
#df_ligands<-left_join(df_ligands,df_resname_change,by=c("ligand"="ligand_name"))
write.csv(df_active_center,"MD_globular_protein_ligand/start/complex.csv",row.names = F)
