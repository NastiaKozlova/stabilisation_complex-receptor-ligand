#prepare files for MD based on chousen structures 
#predicted structures
#and sequense predict structure of protein using robetta server structure_prediction
part_name <- commandArgs(trailingOnly=TRUE)
part_TEMP<-strsplit(part_name,split = ",",fixed = T)[[1]]
part_start<-part_TEMP[1]
protein_name<-part_TEMP[2]
setwd(part_start)
print(protein_name)
if(!dir.exists(paste0(part_start,protein_name))){dir.create(paste0(part_start,protein_name))}
if(!dir.exists(paste0(part_start,protein_name,"/MD_globular_protein/"))){dir.create(paste0(part_start,protein_name,"/MD_globular_protein/"))}
if(!dir.exists(paste0(part_start,protein_name,"/MD_globular_protein/programs/"))){dir.create(paste0(part_start,protein_name,"/MD_globular_protein/programs/"))}
system(command = paste0("cp -r ",part_start,"programs/NAMD_2.14_Linux-x86_64-multicore.tar.gz ",part_start,protein_name,"/MD_globular_protein/programs/"),ignore.stdout=T,wait = T)

system(command = paste0("cp -r ",part_start,"programs/la1.0 ",part_start,protein_name,"/MD_globular_protein/programs/"),ignore.stdout=T,wait = T)
system(command = paste0("cp -r ",part_start,"programs/orient ",part_start,protein_name,"/MD_globular_protein/programs/"),ignore.stdout=T,wait = T)

system(command = paste0("tar -xvzf ",part_start,protein_name,"/MD_globular_protein/programs/NAMD_2.14_Linux-x86_64-multicore.tar.gz -C " ,part_start,protein_name,"/MD_globular_protein/programs/"),ignore.stdout=T,wait = T)

system(command = paste0("rm ",part_start,"r_scripts/MD_globular_protein/start/structure/*"),ignore.stdout=T,wait = T)  
system(command = paste0("cp ",part_start,"start/MD_stabilisation/",protein_name,".pdb ",part_start,"r_scripts/MD_globular_protein/start/structure/",protein_name,".pdb"),ignore.stdout=T,wait = T)
system(command = paste0("cp -r ",part_start,"r_scripts/MD_globular_protein/ ",part_start,protein_name,"/"),ignore.stdout=T,wait = T)
system(command = paste0("cp -r ",part_start,"start/toppar/ ",part_start,protein_name,"/MD_globular_protein/start/"),ignore.stdout=T,wait = T)
