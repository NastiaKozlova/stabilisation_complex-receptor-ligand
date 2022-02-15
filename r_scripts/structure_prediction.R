part_name <- commandArgs(trailingOnly=TRUE)
part_TEMP<-strsplit(part_name,split = ",",fixed = T)[[1]]
part_start<-part_TEMP[1]
v_list_proteins<-part_TEMP[2]
setwd(part_start)
#structure_prediction

if (!dir.exists(paste0(part_start,v_list_proteins))){dir.create(paste0(part_start,v_list_proteins))}
if (!dir.exists(paste0(part_start,v_list_proteins,"/structure_prediction/"))){dir.create(paste0(part_start,v_list_proteins,"/structure_prediction/"))}
if (!dir.exists(paste0(part_start,v_list_proteins,"/structure_prediction/r_scripts/"))){dir.create(paste0(part_start,v_list_proteins,"/structure_prediction/r_scripts/"))}
system(command = paste0("cp -r ",part_start,"start/structure_prediction/pdb/",v_list_proteins,"/ ",part_start,v_list_proteins,"/structure_prediction/pdb/"),ignore.stdout=T,wait = T)
system(command = paste0("cp ",part_start,"start/structure_prediction/sequence/",v_list_proteins,".fasta ",part_start,v_list_proteins,"/structure_prediction/seqs.fasta"),ignore.stdout=T,wait = T)
system(command = paste0("cp ",part_start,"start/sequence/",v_list_proteins,".fasta ",part_start,v_list_proteins,"/structure_prediction/protein_sequence.fasta"),ignore.stdout=T,wait = T)
if(!dir.exists("predicted/")){dir.create("predicted/")}
if(!dir.exists(paste0("predicted/",v_list_proteins))){dir.create(paste0("predicted/",v_list_proteins))}
#system(command = paste0("cp -r ",part_start,"r_scripts/structure_prediction/ ",part_start,v_list_proteins,"/structure_prediction/r_scripts/ "),ignore.stdout=T,wait = T)

Sys.sleep(5)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/structure_prediction/chouse_structure.R ",part_start,v_list_proteins,"/structure_prediction/"),ignore.stdout=T,wait = T)

system(command = paste0("cp -r ",part_start,v_list_proteins,"/structure_prediction/pdb_fin/ ",part_start,"predicted/",v_list_proteins,"/"),ignore.stdout=T,wait = T)
