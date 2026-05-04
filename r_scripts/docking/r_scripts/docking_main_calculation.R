part_protein <- commandArgs(trailingOnly=TRUE)
part_TEMP<-strsplit(part_protein,split = ",",fixed = T)[[1]]
part_start<-part_TEMP[1]
v_list_protein<-part_TEMP[2]
v_center<-part_TEMP[3]
#check protein surface

part_analysis<-paste0(part_start,v_list_protein,"/docking/docking_first/",v_center,"/")
#part_analysis<-paste0(part_name,"docking_first/")
part_scriprs<-paste0(part_start,"r_scripts/docking/r_scripts/")
system(command = paste0("cp ",  part_scriprs,"pdbqt_to_pdbs.py ",part_analysis),ignore.stdout=T,wait = T)
#pdbqt to pdb
system(command = paste0("Rscript --vanilla  ",part_scriprs,"first_docking_start_analysis.R ",part_analysis),ignore.stdout=T,wait = T)
#log to csv
system(command = paste0("Rscript --vanilla  ",part_scriprs,"prepare_log_csv.R ",part_analysis),ignore.stdout=T,wait = T)
#convert all data to the appropriate format for analysis
system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_pre_analysis.R ",part_analysis),ignore.stdout=T,wait = T)
#calculate interactions between receptor and ligands
system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_interactions.R ",part_analysis),ignore.stdout=T,wait = T)

