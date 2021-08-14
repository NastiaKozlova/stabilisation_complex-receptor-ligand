part_start = commandArgs(trailingOnly=TRUE)
setwd(part_start)
part<-paste0(part_start,"structure_prediction/")
if(file.exists(paste0(part,"r_scripts/structure_prediction/chouse_structure.R"))){
  Sys.sleep(5)
}
system(command = paste0("Rscript --vanilla  ",part,"r_scripts/structure_prediction/chouse_structure.R ",part),ignore.stdout=T,wait = T)
Sys.sleep(5)