part_protein <- commandArgs(trailingOnly=TRUE)
part_TEMP<-strsplit(part_protein,split = ",",fixed = T)[[1]]
part_start<-part_TEMP[1]
v_list_protein<-part_TEMP[2]


part_scriprs<-paste0(part_start,"r_scripts/docking/r_scripts/")
part_first<-paste0(part_start,v_list_protein,"/","docking/docking_first/")
part_name<-paste0(part_start,v_list_protein,"/","docking/docking_2/")

if(!dir.exists(part_name)){dir.create(part_name)}
if(!dir.exists(paste0(part_name,"ligand"))){dir.create(paste0(part_name,"ligand"))}
system(command = paste0("cp -r ",part_start,"start/docking/docking_2/ligand_start/ ",part_name),ignore.stdout=T,wait = T)
#ligands_prepare
v_ligand<-list.files(paste0("ligand_start/"))
a<-c()
for (i in 1:length(v_ligand)) {
  b<-strsplit(v_ligand[i],split = ".",fixed = T)[[1]][1]
  a<-c(a,b)
}
v_ligand<-a
df_ligand<-data.frame(matrix(ncol=2,nrow=length(v_ligand)))
colnames(df_ligand)<-c("ligand","C")
df_ligand$ligand<-v_ligand
for (i in 1:nrow(df_ligand)) {
  system(command = paste0("obabel ",part_name,"ligand_start/",df_ligand$ligand[i], ".pdb -O ",part_name,"ligand/",df_ligand$ligand[i], ".pdbqt"),ignore.stdout=T,wait = T)
}

#prepare_receptor

if (!dir.exists(paste0(part_name,"receptor/"))){dir.create(paste0(part_name,"receptor/"))}
if (!dir.exists(paste0(part_name,"active_center/"))){dir.create(paste0(part_name,"active_center/"))}
if (!dir.exists(paste0(part_name,"receptor_start/"))){dir.create(paste0(part_name,"receptor_start/"))}
#df_all_systems<-read.csv(paste0(part_start,"start/all_systems.csv"),stringsAsFactors = F)

system(command = paste0("cp -r ",part_first,"din/complex_structure_center/* ",part_name,"receptor_start/"),ignore.stdout=T,wait = T)
v_receptor<-list.files(paste0(part_name,"receptor_start/"))
df_receptor<-data.frame(matrix(ncol = 2,nrow = length(v_receptor)))
colnames(df_receptor)<-c("receptor","ligand_first")
o<-1
for (o in 1:length(v_receptor)) {
  df_receptor$receptor[o]<-strsplit(v_receptor[o],split = ".",fixed = T)[[1]][1]
  df_receptor$ligand_first[o]<-strsplit(df_receptor$receptor[o],split="_",fixed = T)[[1]][2]
}

i<-1
for (i in 1:nrow(df_receptor)) {
  system(command = paste0(part_start,"programs/MGLTools-1.5.7/bin/pythonsh ",
                          part_start,"programs/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py",
                          " -r ", part_first,"din/complex_structure_center/",df_receptor$receptor[i],".pdb ",
                          " -o ", part_name,"receptor/",df_receptor$receptor[i],".pdbqt ",
                          " -U nphs "))
}
df_cascade_docking<-left_join(df_cascade_docking,df_receptor,by="ligand_first")
df_cascade_docking<-df_cascade_docking%>%mutate(C=NA)
df_active_center<-read.csv(paste0(part_start,"start/active_center.csv"),stringsAsFactors = F)
write.csv(df_active_center,paste0(part_name,"active_center.csv"),row.names=F)
df_active_center<-df_active_center%>%mutate(C=NA)

df_ligand_center<-left_join(df_active_center,df_cascade_docking,by="C")
df_ligand_center<-df_ligand_center%>%mutate(type=center)
df_ligand_center<-df_ligand_center%>%select(type,receptor,ligand)
df_ligand_center<-unique(df_ligand_center)
write.csv(df_ligand_center,paste0(part_name,"ligand_center.csv"),row.names =  F)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_script_2.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/convert_pdbqt_to_pdb.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("chmod +x ",part_name,"convert_pdbqt_to_pdb.py "),ignore.stdout=T,wait = T)
system(command = paste0("python3 ", part_name,"convert_pdbqt_to_pdb.py"),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/convert_log_to_csv.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("chmod +x ",part_name,"prepare_log_csv.py "),ignore.stdout=T,wait = T)
system(command = paste0("python3 ", part_name,"prepare_log_csv.py"),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_pre_analysis.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_group_structure.R ",part_name),ignore.stdout=T,wait = T)
#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/merge_docking_parts.R ",part_name),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_interactions_start.R ",part_name),ignore.stdout=T,wait = T)
#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/analysis.R ",part_name),ignore.stdout=T,wait = T)
