part_prepare = commandArgs(trailingOnly=TRUE)
setwd(part_prepare)
v_ligands<-list.files(paste0("MD_globular_protein_ligand/start/ligands"))
v_receptors<-list.files(paste0("MD_globular_protein_ligand/start/receptor"))
df_active_center<-read.csv(paste0("docking/docking_first/df_all.csv"),stringsAsFactors =  F)
df_active_center<-df_active_center%>%mutate(name=paste0(receptor,"_",ligand,"_",center))

df_ligands<-data.frame(matrix(ncol = 3,nrow = length(v_ligands)))
colnames(df_ligands)<-c("structure_name","complex_name","resname_change")
df_ligands$structure_name<-v_ligands
i<-1
for (i in 1:nrow(df_ligands)) {
  a<-strsplit(df_ligands$structure_name[i],split = "_")[[1]]
  df_ligands$complex_name[i]<-paste0(a[1:(length(a)-1)],collapse =  "_")
}
df_ligands<-left_join(df_ligands,df_active_center,by = c("complex_name"="name"))
i<-1
for (i in 1:nrow(df_ligands)) {
  pdb<-read.pdb(paste0("MD_globular_protein_ligand/start/ligands/",df_ligands$structure_name[i]))
  if(length(unique(pdb$atom$resid))==1){
    df_ligands$resname_change<-unique(pdb$atom$resid)
  }
}
df_ligands<-df_ligands%>%filter(!is.na(resname_change))
write.csv(df_ligands,"MD_globular_protein_ligand/start/complex.csv",row.names = F)
