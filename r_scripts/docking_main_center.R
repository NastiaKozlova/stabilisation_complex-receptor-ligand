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
system(command = paste0("Rscript --vanilla  ",part_scriprs,"prepare_first_docking_center.R ",part_name),ignore.stdout=T,wait = T)
#pdbqt to pdb
system(command = paste0("Rscript --vanilla  ",part_scriprs,"first_docking_start_analysis.R ",part_name),ignore.stdout=T,wait = T)
#log to csv
system(command = paste0("Rscript --vanilla  ",part_scriprs,"prepare_log_csv.R ",part_analysis),ignore.stdout=T,wait = T)
#convert all data to the appropriate format for analysis
system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_pre_analysis_center.R ",part_analysis),ignore.stdout=T,wait = T)
#calculate interactions between receptor and ligands
system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_interactions.R ",part_analysis),ignore.stdout=T,wait = T)
#system(command = paste0("chmod +x ",part_analysis,"prepare_log_csv.py "),ignore.stdout=T,wait = T)
#system(command = paste0("python3 ", part_analysis,"prepare_log_csv.py"),ignore.stdout=T,wait = T)



system(command = paste0("Rscript --vanilla  ",part_scriprs,"RMSD_docking_group_structure.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"calibration_group_structure_center.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_group_center.R ",part_analysis,",2.5"),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_scriprs,"RMSD_merge_docking_center.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"calibration_merge_structure_center.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"merge_docking_center.R ",part_analysis),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_scriprs,"merge_interactions_center.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"affinity_distribution_plot_center.R ",part_analysis),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_scriprs,"complex_structure_center.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"atom_interactions_center.R ",part_analysis),ignore.stdout=T,wait = T)
