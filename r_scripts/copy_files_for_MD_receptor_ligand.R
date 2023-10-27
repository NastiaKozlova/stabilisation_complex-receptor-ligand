part_protein <- commandArgs(trailingOnly=TRUE)
part_TEMP<-strsplit(part_protein,split = ",",fixed = T)[[1]]
part_start<-part_TEMP[1]
v_list_protein<-part_TEMP[2]

part_protein<-paste0(part_start,",",v_list_protein)
part_name<-paste0(part_start,v_list_protein,"/MD_globular_protein_ligand/")
print(v_list_protein)
if(!dir.exists(paste0(part_name))){dir.create(paste0(part_name))}
if(!dir.exists(paste0(part_name,"programs/"))){dir.create(paste0(part_name,"programs/"))}
if(!dir.exists(paste0(part_name,"start/"))){dir.create(paste0(part_name,"start/"))}
if(!dir.exists(paste0(part_name,"start/ligands_surf/"))){dir.create(paste0(part_name,"start/ligands_surf/"))}
if(!dir.exists(paste0(part_name,"start/ligands_center/"))){dir.create(paste0(part_name,"start/ligands_center/"))}
system(command = paste0("cp -r ",part_start,"programs/NAMD_2.14_Linux-x86_64-multicore.tar.gz ",part_name,"programs/"),ignore.stdout=T,wait = T)

if(!dir.exists(paste0(part_name,"start/receptor/"))){dir.create(paste0(part_name,"start/receptor/"))}
system(command = paste0("cp ",part_start,v_list_protein,"/docking/docking_first/din/str_fin_center/* ",part_name,"start/ligands_center/"),ignore.stdout=T,wait = T)

system(command = paste0("cp -r ",part_start,v_list_protein,"/docking/docking_first/receptor_start/start.pdb ", part_name,"start/receptor/start.pdb"))
system(command = paste0("cp -r ",part_start,"programs/NAMD_2.14_Linux-x86_64-multicore.tar.gz ",part_name,"programs/"),ignore.stdout=T,wait = T)
#system(command = paste0("tar -xvzf ",part_name,"programs/NAMD_2.14_Linux-x86_64-multicore.tar.gz ",part_name,"programs/"),ignore.stdout=T,wait = T)


v_ligands<-list.files(paste0(part_start,v_list_protein,"/docking/docking_first/din/structure_merged_center/"))
if (length(v_ligands)>0){
  for (j in 1:length(v_ligands)) {
    system(command = paste0("cp -r ",part_start,v_list_protein,"/docking/docking_first/din/structure_merged_center/",v_ligands[j]," ",part_name,"start/ligands_center/"),ignore.stdout=T,wait = T)
  }
}
v_ligands<-list.files(paste0(part_start,v_list_protein,"/docking/docking_first/din/structure_merged/"))
if (length(v_ligands)>0){
  for (j in 1:length(v_ligands)) {
    system(command = paste0("cp -r ",part_start,v_list_protein,"/docking/docking_first/din/structure_merged/",v_ligands[j]," ",part_name,"start/ligands_surf/"),ignore.stdout=T,wait = T)
  }
}

system(command = paste0("cp -r ",part_start,"start/toppar/ ",part_name,"start/toppar/"),ignore.stdout=T,wait = T)
