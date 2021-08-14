part_start = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
#part_start<-part_analysis
part_name<-part_start
setwd(part_start)
name<-strsplit(part_start,split = "/",fixed = T)[[1]]
based_name<-name[length(name)]
part<-paste0(name[1:(length(name)-1)],collapse = "/")
part_name<-paste0(part,"/")
num_model<-3
#part_start<-paste0(part,"/",based_name[name],'MD/stabilisation/')
name<-1
for (name in 1:length(based_name)) {
  
  setwd(paste0(part_start,"MD/stabilisation/"))
  #create additional directories
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din'))) {dir.create(paste0(part_start,'MD/stabilisation/din'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/tcl'))) {dir.create(paste0(part_start,'MD/stabilisation/din/tcl'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/pdb_second'))) {dir.create(paste0(part_start,'MD/stabilisation/din/pdb_second'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/pdb_second/',based_name[name]))) {dir.create(paste0(part_start,'MD/stabilisation/din/pdb_second/',based_name[name]))}
 
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/pdb_second_complex'))) {dir.create(paste0(part_start,'MD/stabilisation/din/pdb_second_complex'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/pdb_second_complex/',based_name[name]))) {dir.create(paste0(part_start,'MD/stabilisation/din/pdb_second_complex/',based_name[name]))}
  
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/Energy'))) {dir.create(paste0(part_start,'MD/stabilisation/din/Energy'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/SASA'))) {dir.create(paste0(part_start,'MD/stabilisation/din/SASA'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/RMSD'))) {dir.create(paste0(part_start,'MD/stabilisation/din/RMSD'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/RMSF'))) {dir.create(paste0(part_start,'MD/stabilisation/din/RMSF'))}
  #combine MD simulation dcd files
  df_tcl<-data.frame(matrix(nrow = 2,ncol = 2))
  df_tcl[1,1]<-paste0('cd ', part_start,'MD/stabilisation/\npackage require animate')
  
  df_tcl[1,2]<-paste0('mol new {protein/ionized_',based_name[name],'.psf} type {psf}')
  for (i in 1:num_model) {
    if (file.exists(paste0(part_start,'MD/stabilisation/quench/quench_',based_name[name],'_',i,'.dcd'))){
      df_tcl[i+1,1]<-paste0('mol addfile {quench/quench_',based_name[name],'_',i,'.dcd} type {dcd} first 0 last -1 step 10 waitfor all')
    }
  }
  df_tcl[i+2,1]<-paste0('animate write dcd quench/quench_',based_name[name],'.dcd waitfor all')

  df_tcl[i+3,1]<-paste0('mol delete all\n\nexit now')
  
  write.table(df_tcl,file =paste0('din/tcl/combine_',based_name[name],'.tcl'),sep = '\n',na = '' ,row.names = F,col.names = F,quote = F)
  print( paste0('vmd -dispdev text -e ',part_start,'MD/stabilisation/din/tcl/combine_',based_name[name],'.tcl'))
  system(command = paste0('vmd -dispdev text -e ',part_start,'MD/stabilisation/din/tcl/combine_',based_name[name],'.tcl'),ignore.stdout=T,wait = T) 
  
  #RMSD conting
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 13))
  df_tcl[1,1]<-paste0('cd ', part_start,'MD/stabilisation/\n')

  df_tcl[2,1]<-paste0('mol new {protein/ionized_',based_name[name],'.psf} type {psf}')
  df_tcl[2,2]<-paste0('mol addfile {protein/ionized_',based_name[name],'.pdb} type {pdb}') 
  df_tcl[2,3]<-paste0('set protein_0 [atomselect top "protein and name CA"]')
  df_tcl[2,4]<-paste0('mol new {protein/ionized_',based_name[name],'.psf} type {psf}')
  df_tcl[2,5]<-paste0('mol addfile {quench/quench_',based_name[name],'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
  df_tcl[2,6]<-paste0('set n [molinfo top get numframes]')
  df_tcl[2,7]<-paste0('set output [open din/RMSD/',based_name[name],'.txt w] ')
  df_tcl[2,8]<-paste0('for {set i 0} {$i < $n} {incr i} {')
  df_tcl[2,9]<-paste0('set protein_i [atomselect top "protein and name CA" frame $i]')
  df_tcl[2,10]<-paste0('set rmsd [measure rmsd $protein_0 $protein_i]')
  df_tcl[2,11]<-paste0('puts $output "$i $rmsd"')
  df_tcl[2,12]<-paste0('}')
  df_tcl[2,13]<-paste0('puts "output file: $n din/RMSD/',based_name[name],'.txt"')
  df_tcl[2,14]<-paste0('close $output')
  df_tcl[2,15]<-paste0('mol delete all\n\nexit now')
  
  write.table(df_tcl,file =paste0('din/tcl/RMSD_',based_name[name],'.tcl'),sep = '\n',na = '' ,row.names = F,col.names = F,quote = F)
  print( paste0('vmd -dispdev text -e ',part_start,'MD/stabilisation/din/tcl/RMSD_',based_name[name],'.tcl'))
  system(command = paste0('vmd -dispdev text -e ',part_start,'MD/stabilisation/din/tcl/RMSD_',based_name[name],'.tcl'),ignore.stdout=T,wait = T) 
  
  #RMSF counting
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 13))
  
  df_tcl[1,1]<-paste0('cd ', part_start,'MD/stabilisation/')
  df_tcl[2,1]<-paste0('mol new {protein/ionized_',based_name[name],'.psf} type {psf}')
  df_tcl[2,2]<-paste0('mol addfile {quench/quench_',based_name[name],'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
  df_tcl[2,3]<-paste0('set protein [atomselect top "protein and name CA"]')
  df_tcl[2,4]<-paste0('set n [molinfo top get numframes]')
  df_tcl[2,5]<-paste0('set output [open din/RMSF/',based_name[name],'.txt w] ')
  df_tcl[2,6]<-paste0('set rmsf [measure rmsf $protein]')
  df_tcl[2,7]<-paste0('foreach x $rmsf {')
  df_tcl[2,8]<-paste0('puts $output $x')
  df_tcl[2,9]<-paste0('}')
  df_tcl[2,10]<-paste0('puts "output file: $n din/RMSF/',based_name[name],'.txt"')
  df_tcl[2,11]<-paste0('close $output')
  df_tcl[2,12]<-paste0('mol delete all\n\nexit now')
  write.table(df_tcl,file =paste0('din/tcl/RMSF_',based_name[name],'.tcl'),sep = '\n',na = '' ,row.names = F,col.names = F,quote = F)
  system(command = paste0('vmd -dispdev text -e ',part_start,'MD/stabilisation/din/tcl/RMSF_',based_name[name],'.tcl'),ignore.stdout=T,wait = T) 
  
  #Secondary structure counting
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 7))
  df_tcl[1,1]<-paste0('cd ', part_start,'MD/stabilisation/')
  df_tcl[2,1]<-paste0('mol new {protein/ionized_',based_name[name],'.psf} type {psf}')
  df_tcl[2,2]<-paste0('mol addfile {quench/quench_',based_name[name],'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
  df_tcl[2,3]<-paste0('set nf [molinfo top get numframes]')
  df_tcl[2,4]<-paste0('for {set i 0 } {$i < $nf} {incr i} {')
  df_tcl[2,5]<-paste0('[atomselect top "protein" frame $i] writepdb din/pdb_second/',based_name[name],'/frame_$i.pdb')
  df_tcl[2,6]<-paste0('}')
  df_tcl[2,7]<-'mol delete all\n\nexit now'
  write.table(df_tcl,file =paste0('din/tcl/Second_str_',based_name[name],'.tcl'),sep = '\n',na = '' ,row.names = F,col.names = F,quote = F)
  print(paste0('vmd -dispdev text -e ',part_start,'MD/stabilisation/din/tcl/Second_str_',based_name[name],'.tcl'))
  system(command = paste0('vmd -dispdev text -e ',part_start,'MD/stabilisation/din/tcl/Second_str_',based_name[name],'.tcl'),ignore.stdout=T,wait = T) 
  
  #Secondary structure counting
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 7))
  df_tcl[1,1]<-paste0('cd ', part_start,'MD/stabilisation/')
  df_tcl[2,1]<-paste0('mol new {protein/ionized_',based_name[name],'.psf} type {psf}')
  df_tcl[2,2]<-paste0('mol addfile {quench/quench_',based_name[name],'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
  df_tcl[2,3]<-paste0('set nf [molinfo top get numframes]')
  df_tcl[2,4]<-paste0('for {set i 0 } {$i < $nf} {incr i} {')
  df_tcl[2,5]<-paste0('[atomselect top "not (water or ions)" frame $i] writepdb din/pdb_second_complex/',based_name[name],'/frame_$i.pdb')
  df_tcl[2,6]<-paste0('}')
  df_tcl[2,7]<-'mol delete all\n\nexit now'
  write.table(df_tcl,file =paste0('din/tcl/Second_str_',based_name[name],'_complex.tcl'),sep = '\n',na = '' ,row.names = F,col.names = F,quote = F)
  print(paste0('vmd -dispdev text -e ',part_start,'MD/stabilisation/din/tcl/Second_str_',based_name[name],'.tcl'))
  system(command = paste0('vmd -dispdev text -e ',part_start,'MD/stabilisation/din/tcl/Second_str_',based_name[name],'_complex.tcl'),ignore.stdout=T,wait = T) 
  

  #Energy counting
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  df_tcl[1,1]<-paste0('cd ', part_start,'MD/stabilisation/\npackage require namdenergy')
  df_tcl[1,2]<-paste0('mol new {protein/ionized_',based_name[name],'.psf} type {psf}')
  df_tcl[1,3]<-paste0('mol addfile {quench/quench_',based_name[name],'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
  df_tcl[1,4]<-paste0('set sel2 [atomselect top "protein"]')
  df_tcl[1,5]<-paste0('set sel1 [atomselect top "not (protein or water or ions)"]')
  df_tcl[1,6]<-paste0('\nnamdenergy -sel $sel1 $sel2 -vdw -elec -nonb -cutoff 12 -skip 0 -ofile din/Energy/interactions_',based_name[name],'.txt -switch 10 -exe /home/nastia/NAMD_2.14_Linux-x86_64-multicore/namd2 -par ',part_name,'start/toppar/par_all36_carb.prm -par ',part_name,'start/toppar/par_all36_cgenff.prm -par ',part_name,'start/toppar/par_all36_lipid.prm -par ',part_name,'start/toppar/par_all36m_prot.prm -par ',part_name,'start/toppar/par_all36_na.prm -par ',part_name,'start/toppar/par_all36_prot.prm -par ',part_name,'start/toppar/toppar_water_ions_namd.str')
  
  df_tcl[1,7]<-paste0('\nnamdenergy -sel $sel2  -bond -angl -dihe -impr -conf -vdw -elec -nonb -all -cutoff 12 -skip 0 -ofile din/Energy/protein_',based_name[name],'.txt -switch 10 -exe /home/nastia/NAMD_2.14_Linux-x86_64-multicore/namd2 -par ',part_name,'start/toppar/par_all36_carb.prm -par ',part_name,'start/toppar/par_all36_cgenff.prm -par ',part_name,'start/toppar/par_all36_lipid.prm -par ',part_name,'start/toppar/par_all36m_prot.prm -par ',part_name,'start/toppar/par_all36_na.prm -par ',part_name,'start/toppar/par_all36_prot.prm -par ',part_name,'start/toppar/toppar_water_ions_namd.str')
  df_tcl[1,8]<-'mol delete all'
  df_tcl[1,9]<-'\n\nexit now'
  write.table(df_tcl,file =paste0('din/tcl/Energy_',based_name[name],'.tcl'),sep = '\n',na = '' ,row.names = F,col.names = F,quote = F)
  print(paste0('vmd -dispdev text -e ',part_start,'MD/stabilisation/din/tcl/Energy_',based_name[name],'.tcl'))
  system(command = paste0('vmd -dispdev text -e ',part_start,'MD/stabilisation/din/tcl/Energy_',based_name[name],'.tcl'),ignore.stdout=T,wait = T) 
  
  #SASA counting
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 9))
  df_tcl[1,1]<-paste0('cd ', part_start,'MD/stabilisation/\npackage require measure')
  df_tcl[1,2]<-paste0('mol new {protein/ionized_',based_name[name],'.psf} type {psf}')
  df_tcl[1,3]<-paste0('mol addfile {quench/quench_',based_name[name],'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
  df_tcl[1,4]<-paste0('set protein [atomselect top "protein"]')
  df_tcl[1,5]<-paste0('set n [molinfo top get numframes]')
  df_tcl[1,6]<-paste0('set output [open din/SASA/',based_name[name],'.txt w] ')
  df_tcl[1,7]<-paste0('for {set i 0} {$i < $n} {incr i} {')
  df_tcl[1,8]<-paste0('molinfo top set frame $i')
  df_tcl[1,9]<-paste0('set sasa_protein [measure sasa 1.4 $protein]')
  df_tcl[1,10]<-paste0('puts $output " $i $sasa_protein "')
  df_tcl[1,11]<-paste0('}')
  df_tcl[1,12]<-paste0('puts "output file: $n din/SASA/',based_name[name],'.txt"')
  df_tcl[1,13]<-paste0('close $output')
  df_tcl[1,14]<-paste0('mol delete all\n\nexit now')
  write.table(df_tcl,file =paste0('din/tcl/SASA_',based_name[name],'.tcl'),sep = '\n',na = '' ,row.names = F,col.names = F,quote = F)
  print(paste0('vmd -dispdev text -e ',part_start,'MD/stabilisation/din/tcl/SASA_',based_name[name],'.tcl'))
  system(command = paste0('vmd -dispdev text -e ',part_start,'MD/stabilisation/din/tcl/SASA_',based_name[name],'.tcl'),ignore.stdout=T,wait = T) 
}

