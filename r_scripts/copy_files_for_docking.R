part_name <- commandArgs(trailingOnly=TRUE)
part_TEMP<-strsplit(part_name,split = ",",fixed = T)[[1]]
part_start<-part_TEMP[1]
v_list_protein<-part_TEMP[2]
part_name<-  paste0(part_start,v_list_protein,"/docking/")
setwd(part_start)
#copying sctipts for docking
system(command = paste0("cp -r ",part_start,"r_scripts/docking/ ",part_start,v_list_protein,"/"),ignore.stdout=T,wait = T)
#create directories
if(!dir.exists(paste0(part_name))){dir.create(paste0(part_name))}
if(!dir.exists(paste0(part_name,"docking_first/"))){dir.create(paste0(part_name,"docking_first/"))}
if(!dir.exists(paste0(part_name,"docking_first/center"))){dir.create(paste0(part_name,"docking_first/center"))}
if(!dir.exists(paste0(part_name,"docking_first/surf"))){dir.create(paste0(part_name,"docking_first/surf"))}
if(!dir.exists(paste0(part_name,"docking_first/center/receptor_start/"))){dir.create(paste0(part_name,"docking_first/center/receptor_start/"))}
if(!dir.exists(paste0(part_name,"docking_first/surf/receptor_start/"))){dir.create(paste0(part_name,"docking_first/surf/receptor_start/"))}

if(!dir.exists(paste0(part_name,"docking_first/center/receptor/"))){dir.create(paste0(part_name,"docking_first/center/receptor/"))}
if(!dir.exists(paste0(part_name,"docking_first/surf/receptor/"))){dir.create(paste0(part_name,"docking_first/surf/receptor/"))}
#copying ligand's structures for docking
system(command = paste0("cp -r ",part_start,"start/docking/docking_first/ligand/ ",part_name,"/docking_first/center/"),ignore.stdout=T,wait = T)
system(command = paste0("cp -r ",part_start,"start/docking/docking_first/ligand/ ",part_name,"/docking_first/surf/"),ignore.stdout=T,wait = T)
#copying ligand field search 
system(command = paste0("cp -r ",part_start,"start/docking/ligand_field.csv ",part_name,"/docking_first/center/"),ignore.stdout=T,wait = T)
system(command = paste0("cp -r ",part_start,"start/docking/ligand_field.csv ",part_name,"/docking_first/surf/"),ignore.stdout=T,wait = T)

#copying receptor structure for docking
system(command = paste0("cp -r ",part_start,"output/stabilisation/",v_list_protein,".pdb ",part_name,"docking_first/center/receptor_start/start.pdb"),ignore.stdout=T,wait = T)
system(command = paste0("cp -r ",part_start,"output/stabilisation/",v_list_protein,".pdb ",part_name,"docking_first/surf/receptor_start/start.pdb"),ignore.stdout=T,wait = T)

#copying amino acids for active center for docking
system(command = paste0("cp ",part_start,"start/docking/active_center/",v_list_protein,".csv ",part_name,"/docking_first/center/active_center.csv"),ignore.stdout=T,wait = T)
#copying programs for docking
system(command = paste0("cp -r ",part_start,"programs/ ",part_name,"docking_first/center/"),ignore.stdout=T,wait = T)
system(command = paste0("cp -r ",part_start,"programs/ ",part_name,"docking_first/surf/"),ignore.stdout=T,wait = T)
