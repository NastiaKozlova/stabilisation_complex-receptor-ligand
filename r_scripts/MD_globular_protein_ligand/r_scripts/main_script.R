#name your partitions part and part_start 
#put your pdb files into start/structure
part<-"part to the your directory"
part<-"/home/nastia/mem/MD_globular_protein_ligand/"

part_start<-part

library(dplyr)
library(bio3d)

#change v_namd in the script to the your namd path
#structure_name in "start/complex.csv"shoud not contain .pdb
#does not support modified proteins
#to run modify protein run first prepare_to_stabilisation_MD.R and swap pdb and psf structure of modify protein chains in 
#prepared_structures
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_to_stabilisation_MD.R ",part_start),ignore.stdout=T,wait = T)
#system(command = paste0("vmd -dispdev text -e ",part_start,"r_scripts/prepare_MD.tcl "),ignore.stdout=T,wait = T) 
#namd
print(paste0("#run namd2 from "))
print(paste0(part_start,"r_scripts/namd_script.txt"))

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_tcl_din.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/second_stucture_compare.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/Ramachadran.R ",part_start),ignore.stdout=T,wait = T)
#make plot
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/make_plots_RMSD_RMSF.R ",part_start),ignore.stdout=T,wait = T)
