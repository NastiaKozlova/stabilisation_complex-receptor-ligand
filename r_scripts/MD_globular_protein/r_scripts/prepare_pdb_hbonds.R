part_name = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
#part_name<-part_name

setwd(part_name)
name<-list.files('start/structure/')

a<-c()
for(i in 1:length(name)){
  b<-strsplit(name[i],split = '.',fixed = T)[[1]][1]
  a<-c(a,b)
}
based_name<-a
name<-1
part<-part_name
num_model<-100

for (name in 1:length(based_name)) {
  part<-paste0(part_name,based_name[name],'/MD/stabilisation/')
  setwd(paste0(part))
  #create additional directories
#  if (!dir.exists(paste0(part,'din'))) {dir.create(paste0(part,'din'))}
#  if (!dir.exists(paste0(part,'din/tcl'))) {dir.create(paste0(part,'din/tcl'))}
  if (!dir.exists(paste0(part,'din/pdb_second'))) {dir.create(paste0(part,'din/pdb_second'))}
#  if (!dir.exists(paste0(part,'din/pdb_second/',based_name[name]))) {dir.create(paste0(part,'din/pdb_second/',based_name[name]))}
  
  if (!dir.exists(paste0(part,'din/pdb_second/hbonds_',based_name[name]))) {dir.create(paste0(part,'din/pdb_second/hbonds_',based_name[name]))}
  
#  if (!dir.exists(paste0(part,'din/Energy'))) {dir.create(paste0(part,'din/Energy'))}
#  if (!dir.exists(paste0(part,'din/SASA'))) {dir.create(paste0(part,'din/SASA'))}
#  if (!dir.exists(paste0(part,'din/RMSD'))) {dir.create(paste0(part,'din/RMSD'))}
#  if (!dir.exists(paste0(part,'din/RMSF'))) {dir.create(paste0(part,'din/RMSF'))}
  #Secondary structure counting
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 7))
  df_tcl[1,1]<-paste0('cd ', part,'')
  df_tcl[2,1]<-paste0('mol new {protein/ionized_',based_name[name],'.psf} type {psf}')
  df_tcl[2,2]<-paste0('mol addfile {quench/quench_',based_name[name],'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
  df_tcl[2,3]<-paste0('set nf [molinfo top get numframes]')
  df_tcl[2,4]<-paste0('for {set i 0 } {$i < $nf} {incr i} {')
  df_tcl[2,5]<-paste0('[atomselect top all frame $i] writepdb din/pdb_second/hbonds_',based_name[name],'/frame_$i.pdb')
  df_tcl[2,6]<-paste0('}')
  df_tcl[2,7]<-'mol delete all\n\nexit now'
  write.table(df_tcl,file =paste0('din/tcl/Second_str_hbonds_',based_name[name],'.tcl'),sep = '\n',na = '' ,row.names = F,col.names = F,quote = F)
  print(paste0('vmd -dispdev text -e ',part,'din/tcl/Second_str_hbonds_',based_name[name],'.tcl'))
  system(command = paste0('vmd -dispdev text -e ',part,'din/tcl/Second_str_hbonds_',based_name[name],'.tcl'),ignore.stdout=T,wait = T) 
}

