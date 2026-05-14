part_protein <- commandArgs(trailingOnly=TRUE)
part_TEMP<-strsplit(part_protein,split = ",",fixed = T)[[1]]
part_start<-part_TEMP[1]
v_list_protein<-part_TEMP[2]
v_center<-part_TEMP[3]

part_analysis<-paste0(part_start,v_list_protein,"/docking/docking_first/",v_center,"/")

part_scriprs<-paste0(part_start,"r_scripts/docking/r_scripts/")

system(command = paste0("Rscript --vanilla  ",part_scriprs,"RMSD_merge_docking.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"calibration_merge_structure_center.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"merge_docking_structures.R ",part_analysis,",",7),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_scriprs,"merge_interactions_center.R ",part_analysis),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_scriprs,"complex_structure_center.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"atom_interactions_center.R ",part_analysis),ignore.stdout=T,wait = T)