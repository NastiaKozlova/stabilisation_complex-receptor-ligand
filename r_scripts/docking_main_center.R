part_protein <- commandArgs(trailingOnly=TRUE)
part_TEMP<-strsplit(part_protein,split = ",",fixed = T)[[1]]
part_start<-part_TEMP[1]
v_list_protein<-part_TEMP[2]
#check protein surface

part_name<-paste0(part_start,v_list_protein,"/docking/")
part_analysis<-paste0(part_name,"docking_first/")
part_scriprs<-paste0(part_start,"r_scripts/docking/r_scripts/")

system(command = paste0("cp -r ",part_start,"start/toppar/ ",part_analysis),ignore.stdout=T,wait = T)
#copying amino acids for active center for docking
system(command = paste0("cp ",part_start,"start/docking/active_center/",v_list_protein,".csv ",part_name,"/docking_first/active_center.csv"),ignore.stdout=T,wait = T)

#run docking
system(command = paste0("Rscript --vanilla  ",part_scriprs,"prepare_first_docking_main.R ",part_name),ignore.stdout=T,wait = T)
#pdbqt to pdb
system(command = paste0("Rscript --vanilla  ",part_scriprs,"first_docking_start_analysis.R ",part_name),ignore.stdout=T,wait = T)
#log to csv
system(command = paste0("Rscript --vanilla  ",part_scriprs,"prepare_log_csv.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("chmod +x ",part_analysis,"prepare_log_csv.py "),ignore.stdout=T,wait = T)
system(command = paste0("python3 ", part_analysis,"prepare_log_csv.py"),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_pre_analysis.R ",part_analysis),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_group_structure.R ",part_analysis),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_interactions_center.R ",part_analysis),ignore.stdout=T,wait = T)
