part_name <- commandArgs(trailingOnly=TRUE)
part_TEMP<-strsplit(part_name,split = ",",fixed = T)[[1]]
part_start<-part_TEMP[1]
v_list_protein<-part_TEMP[2]
part_name<-  paste0(part_start,v_list_protein,"/docking/")
setwd(part_start)
#copying sctipts for docking
system(command = paste0("cp -r ",part_start,"r_scripts/docking/ ",part_start,v_list_protein,"/"),ignore.stdout=T,wait = T)
#create directories
if(!dir.exists(paste0(part_name,"docking_first/"))){dir.create(paste0(part_name,"docking_first/"))}
if(!dir.exists(paste0(part_name,"docking_first/receptor_start/"))){dir.create(paste0(part_name,"docking_first/receptor_start/"))}
if(!dir.exists(paste0(part_name,"docking_first/receptor/"))){dir.create(paste0(part_name,"docking_first/receptor/"))}
#copying ligand's structures for docking
system(command = paste0("cp -r ",part_start,"start/docking/docking_first/ligand/ ",part_name,"/docking_first/"),ignore.stdout=T,wait = T)
#copying receptor structure for docking
system(command = paste0("cp -r ",part_start,"output/stabilisation/fin_structure/",v_list_protein,".pdb ",part_name,"/docking_first/receptor_start/start.pdb"),ignore.stdout=T,wait = T)
#copying amino acids for active center for docking
system(command = paste0("cp ",part_start,"start/docking/active_center/",v_list_protein,"/active_center.csv ",part_name,"/docking_first/active_center.csv"),ignore.stdout=T,wait = T)
#copying programs for docking
system(command = paste0("cp -r ",part_start,"programs/ ",part_name,"docking_first/"),ignore.stdout=T,wait = T)
