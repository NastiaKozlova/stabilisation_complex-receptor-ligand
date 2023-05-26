part_protein <- commandArgs(trailingOnly=TRUE)
part_TEMP<-strsplit(part_protein,split = ",",fixed = T)[[1]]
part_start<-part_TEMP[1]
v_list_protein<-part_TEMP[2]
#check protein surface

part_name<-paste0(part_start,v_list_protein,"/docking/")
part_analysis<-paste0(part_name,"docking_first/")
part_scriprs<-paste0(part_start,"r_scripts/docking/r_scripts/")
  
system(command = paste0("cp -r ",part_start,"start/toppar/ ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"check_surface.R ",part_name),ignore.stdout=T,wait = T)
#add surf active center, optional 
system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_add_serf_active_centers.R ",part_name),ignore.stdout=T,wait = T)
#run docking
system(command = paste0("Rscript --vanilla  ",part_scriprs,"prepare_first_docking_main.R ",part_name),ignore.stdout=T,wait = T)
#pdbqt to pdb
system(command = paste0("Rscript --vanilla  ",part_scriprs,"first_docking_start_analysis.R ",part_name),ignore.stdout=T,wait = T)
#log to csv
system(command = paste0("Rscript --vanilla  ",part_scriprs,"prepare_log_csv.R ",part_analysis),ignore.stdout=T,wait = T)
#convert all data to the appropriate format for analysis
system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_pre_analysis.R ",part_analysis),ignore.stdout=T,wait = T)
#calculate interactions between receptor and ligands
system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_interactions.R ",part_analysis),ignore.stdout=T,wait = T)
#calculate RMSD between all structures in the groups
system(command = paste0("Rscript --vanilla  ",part_scriprs,"RMSD_docking_group_structure.R ",part_analysis),ignore.stdout=T,wait = T)
#make a calibration graph of RMSD cutoff
system(command = paste0("Rscript --vanilla  ",part_scriprs,"calibration_group_structure.R ",part_analysis),ignore.stdout=T,wait = T)
#sort structures in the groups
system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_group_structure.R ",part_analysis,",",1),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_scriprs,"RMSD_merge_docking_parts.R ",part_analysis),ignore.stdout=T,wait = T)
#make a calibration graph of RMSD cutoff
system(command = paste0("Rscript --vanilla  ",part_scriprs,"calibration_merge_structure.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"merge_docking_parts.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"merge_interactions.R ",part_analysis),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_scriprs,"complex_structure_surf.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"atom_interactions_surf.R ",part_analysis),ignore.stdout=T,wait = T)
