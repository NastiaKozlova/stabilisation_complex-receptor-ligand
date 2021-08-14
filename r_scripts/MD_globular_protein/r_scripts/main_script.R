#name your partitions part and part_start 
#put your pdb files into start/structure
part_start = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(bio3d)
#change v_namd in the script to the your namd path
#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_to_stabilisation_MD.R ",part_start),ignore.stdout=T,wait = T)
#namd
print(paste0("#run namd2 from "))
print(paste0(part_start,"r_scripts/namd_script.txt"))

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_tcl_din.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/second_stucture_compare.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/Ramachadran.R ",part_start),ignore.stdout=T,wait = T)
#make plot
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/make_plots_RMSD_RMSF.R ",part_start),ignore.stdout=T,wait = T)
