part_protein <- commandArgs(trailingOnly=TRUE)
part_TEMP<-strsplit(part_protein,split = ",",fixed = T)[[1]]
part_start<-part_TEMP[1]
v_list_protein<-part_TEMP[2]
#check protein surface

#part_name<-paste0(part_start,v_list_protein,"/docking/")
part_analysis<-paste0(part_start,v_list_protein,"/docking/docking_first/surf/")
part_scriprs<-paste0(part_start,"r_scripts/docking/r_scripts/")

system(command = paste0("cp -r ",part_start,"start/toppar/ ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"check_surface.R ",part_analysis),ignore.stdout=T,wait = T)
#add surf active center, optional 
system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_add_serf_active_centers.R ",part_analysis),ignore.stdout=T,wait = T)
#run docking
#system(command = paste0("Rscript --vanilla  ",part_scriprs,"prepare_first_docking_main.R ",part_name),ignore.stdout=T,wait = T)