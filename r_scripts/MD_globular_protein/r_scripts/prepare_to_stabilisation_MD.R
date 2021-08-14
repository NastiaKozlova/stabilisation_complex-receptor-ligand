part_start = commandArgs(trailingOnly=TRUE)
# namd run script 
v_namd<-"namd run script"
#quantity of 1 ns MD simulations, change it to Modify length of MD simulation 
num_din<-100
library(dplyr)
library(bio3d)
setwd(part_start)
name<-list.files("start/structure/")

a<-c()
for(i in 1:length(name)){
  b<-strsplit(name[i],split = ".",fixed = T)[[1]][1]
  a<-c(a,b)
}
name<-a
part<-1
for (part in 1:length(name)) {
  print(paste0(name[part]))
  if(!dir.exists(paste0(name[part]))){dir.create(paste0(name[part]))}
  if(!dir.exists(paste0(name[part],"/MD"))){dir.create(paste0(name[part],"/MD"))}
  if(!dir.exists(paste0(name[part],"/MD/stabilisation"))){dir.create(paste0(name[part],"/MD/stabilisation"))}
  if(!dir.exists(paste0(name[part],"/MD/stabilisation/protein"))){dir.create(paste0(name[part],"/MD/stabilisation/protein"))}
  if(!dir.exists(paste0(name[part],"/MD/stabilisation/pdb"))){dir.create(paste0(name[part],"/MD/stabilisation/pdb"))}
  if(!dir.exists(paste0(name[part],"/MD/stabilisation/quench"))){dir.create(paste0(name[part],"/MD/stabilisation/quench"))}
  if(!dir.exists(paste0(name[part],"/MD/stabilisation/dcd"))){dir.create(paste0(name[part],"/MD/stabilisation/dcd"))}
  pdb<-read.pdb(paste0(part_start,"start/structure/",name[part],".pdb"))
  write.pdb(pdb,paste0(part_start,name[part],'/MD/stabilisation/protein/start.pdb'))
  
  #prepare psf and pdb to run NAMD
  
  df_psfgen<-data.frame(matrix(ncol = 1,nrow = 1))
  df_psfgen[1,1]<-paste0('cd ',part_start,name[part],'/MD/stabilisation/protein
  mol delete all
  package require psfgen 
  lappend auto_path ',part_start,'programs/la1.0
  lappend auto_path ',part_start,'programs/orient
  package require Orient
  namespace import Orient::orient
  resetpsf
  topology ',part_start,'start/toppar/top_all36_prot.rtf
  topology ',part_start,'start/toppar/toppar_water_ions_namd.str
  
  pdbalias residue HIS HSE
  pdbalias atom ILE CD1 CD
  segment U { pdb start.pdb
  }',
                         '\ncoordpdb start.pdb U
  regenerate angles dihedrals
  guesscoord
    
  writepdb ',name[part],'_TEMP.pdb
  writepsf ',name[part],'.psf
    
  mol delete all
  mol new ',name[part],'.psf
  mol addfile ',name[part],'_TEMP.pdb
  set sel [atomselect top all]
  set gec [measure center $sel]
  $sel moveby [vecscale -1.0 $gec]
  set I [draw principalaxes $sel]
  set A [orient $sel [lindex $I 2] {0 0 1}]
  $sel move $A
  $sel writepdb ',name[part],'.pdb
  
  file delete ',name[part],'_TEMP.pdb
  mol delete all\n\nexit now')
  write.table(df_psfgen,paste0(name[part],'/MD/stabilisation/protein/psfgen_',name[part],'.tcl'),col.names = F,row.names = F,quote = F)
  system(command = paste0("vmd -dispdev text -e ",part_start,name[part],'/MD/stabilisation/protein/psfgen_',name[part],'.tcl'),ignore.stdout=T,wait = T) 
}
for (part in 1:length(name)) {
  print(paste0(name[part]))
  pdb_start<-read.pdb(paste0(name[part],'/MD/stabilisation/protein/',name[part],'.pdb'))
  df_pdb<-pdb_start$atom
  
  x_min<-min(df_pdb$x)-20
  y_min<-min(df_pdb$y)-20
  z_min<-min(df_pdb$z)-20
  
  x_max<-max(df_pdb$x)+20
  y_max<-max(df_pdb$y)+20
  z_max<-max(df_pdb$z)+20
  #prepare script to solvate and ionize protein pstructure
  df_psfgen<-data.frame(matrix(ncol = 1,nrow = 1))
  df_psfgen[1,1]<-paste0('cd ',part_start,name[part],'/MD/stabilisation/protein\n',
                         'package require solvate \n','package require autoionize \n',
                         'solvate ',name[part],'.psf ',name[part],'.pdb -o solvate_',name[part],' -b 1.5 -minmax {{',x_min,' ',y_min,' ',z_min,'} {',x_max,' ',y_max,' ',z_max,'}}\n',
                         'autoionize -psf solvate_',name[part],'.psf -pdb solvate_',name[part],'.pdb -sc 0.15 -o ionized_',name[part],
                         '\nmol delete all\n\nexit now')
  write.table(df_psfgen,paste0(name[part],'/MD/stabilisation/protein/solvate_',name[part],'.tcl'),col.names = F,row.names = F,quote = F)
  system(command = paste0("vmd -dispdev text -e ",part_start,name[part],'/MD/stabilisation/protein/solvate_',name[part],'.tcl'),ignore.stdout=T,wait = T) 
  
  x_mean<-mean(x_max+x_min)/2
  y_mean<-mean(y_max+y_min)/2
  z_mean<-mean(z_max+z_min)/2
  #energy minimisation 
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  df_tcl[1,1]<-paste0('structure protein/ionized_',name[part],'.psf',
                      '\ncoordinates protein/ionized_',name[part],'.pdb',
                      '\nset temperature 300',
                      '\nfirsttimestep 0',
                      '\n#############################################################',
                      '\n## SIMULATION PARAMETERS ##',
                      '\n#############################################################',
                      '\nparaTypeCharmm on',
                      '\nparameters ',part_start,'start/toppar/par_all36_carb.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_cgenff.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_lipid.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36m_prot.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_na.prm',
                      '\nparameters ',part_start,'start/toppar/toppar_water_ions_namd.str',
                      '\ntemperature $temperature',
                      '\n# Force-Field Parameters',
                      '\nexclude scaled1-4',
                      '\n1-4scaling 1.0',
                      '\nswitching on',
                      '\nswitchdist 8.0',
                      '\ncutoff 12.0',
                      '\npairlistdist 13.5',
                      '\nstepspercycle 20',
                      '\nPME on',
                      '\nPMETolerance 0.000001',
                      '\nPMEGridSizeX ',round(x_max-x_min,digits = 0), 
                      '\nPMEGridSizeY ',round(y_max-y_min,digits = 0),	
                      '\nPMEGridSizeZ ',round(z_max-z_min,digits = 0),
                      '\nfullElectFrequency 4',
                      '\ncellBasisVector1 ',round(x_max-x_min,digits = 0),' 0.0 0.0',
                      '\ncellBasisVector2 0.0 ',round(y_max-y_min,digits = 0),' 0.0',
                      '\ncellBasisVector3 0.0 0.0 ',round(z_max-z_min,digits = 0),
                      '\ncellOrigin ',x_mean,' ',y_mean,' ',z_mean,
                      '\nwrapAll on',
                      '\nminimization on',
                      '\nnumsteps 100000',
                      '\noutputenergies 100',
                      '\ndcdfreq 100',
                      '\ndcdfile dcd/min_',name[part],'.dcd',
                      '\nbinaryoutput no',
                      '\noutputname pdb/min_',name[part])
  write.table(df_tcl,paste0(name[part],'/MD/stabilisation/min_',name[part],'.conf'),row.names = F,col.names = F,quote = F)
  
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  #structure heating 
  df_tcl[1,1]<-paste0('structure protein/ionized_',name[part],'.psf',
                      '\ncoordinates pdb/min_',name[part],'.coor',
                      '\n#############################################################',
                      '\n## SIMULATION PARAMETERS ##',
                      '\n#############################################################',
                      '\nparaTypeCharmm on',
                      '\nparameters ',part_start,'start/toppar/par_all36_prot.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_carb.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_cgenff.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_lipid.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_na.prm',
                      '\nparameters ',part_start,'start/toppar/toppar_water_ions_namd.str',
                      '\n# Force-Field Parameters',
                      '\nexclude scaled1-4',
                      '\n1-4scaling 1.0',
                      '\nswitching on',
                      '\nswitchdist 8.0',
                      '\ncutoff 12.0',
                      '\npairlistdist 13.5',
                      '\nstepspercycle 20',
                      '\nPME on',
                      '\nPMETolerance 0.000001',
                      '\nPMEGridSizeX ',round(x_max-x_min,digits = 0),
                      '\nPMEGridSizeY ',round(y_max-y_min,digits = 0),	
                      '\nPMEGridSizeZ ',round(z_max-z_min,digits = 0),
                      '\nfullElectFrequency 4',
                      '\ncellBasisVector1 ',round(x_max-x_min,digits = 0),' 0.0 0.0',
                      '\ncellBasisVector2 0.0 ',round(y_max-y_min,digits = 0),' 0.0',
                      '\ncellBasisVector3 0.0 0.0 ',round(z_max-z_min,digits = 0),
                      '\ncellOrigin ',x_mean,' ',y_mean,' ',z_mean,
                      '\nwrapAll on',
                      '\ntimestep 1.0',
                      '\ntemperature 0',
                      '\nreassignFreq 10',
                      '\nreassignIncr 0.01',
                      '\nreassignHold 300',
                      '\nfirsttimestep 0',
                      '\nnumsteps 300000',
                      '\noutputenergies 100',
                      '\ndcdfreq 100',
                      '\ndcdfile dcd/head_',name[part],'.dcd',
                      '\nbinaryoutput no',
                      '\noutputname pdb/head_',name[part])
  write.table(df_tcl,paste0(name[part],'/MD/stabilisation/heat_',name[part],'.conf'),row.names = F,col.names = F,quote = F)  
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  #first eqvilibration
  df_tcl[1,1]<-paste0('structure protein/ionized_',name[part],'.psf',
                      '\ncoordinates pdb/head_',name[part],'.coor',
                      '\n#############################################################',
                      '\n## SIMULATION PARAMETERS ##',
                      '\n#############################################################',
                      '\nparaTypeCharmm on',
                      '\nparameters ',part_start,'start/toppar/par_all36_prot.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_carb.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_cgenff.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_lipid.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_na.prm',
                      '\nparameters ',part_start,'start/toppar/toppar_water_ions_namd.str',
                      '\nexclude		scaled1-4',
                      '\n1-4scaling	1.0',
                      '\nswitching	on',
                      '\nswitchdist	8.0',
                      '\ncutoff		12.0',
                      '\npairlistdist	13.5',
                      '\nstepspercycle	20',
                      '\nPME		on',
                      '\nPMETolerance 	0.000001',
                      '\nPMEGridSizeX ',round(x_max-x_min,digits = 0),
                      '\nPMEGridSizeY ',round(y_max-y_min,digits = 0),	
                      '\nPMEGridSizeZ ',round(z_max-z_min,digits = 0),
                      '\nfullElectFrequency 4',
                      '\ncellBasisVector1 ',round(x_max-x_min,digits = 0),' 0.0 0.0',
                      '\ncellBasisVector2 0.0 ',round(y_max-y_min,digits = 0),' 0.0',
                      '\ncellBasisVector3 0.0 0.0 ',round(z_max-z_min,digits = 0),
                      '\nwrapAll on',
                      '\ntimestep		1.0',
                      '\ntemperature		300',
                      '\nreassignFreq		10',
                      '\nnumsteps		2000000',
                      '\nseed			11552514',
                      '\noutputenergies	100',
                      '\ndcdfreq 100',
                      '\ndcdfile dcd/eqv_',name[part],'.dcd',
                      '\nbinaryoutput no',
                      '\noutputname pdb/eqv_',name[part])
  
  write.table(df_tcl,paste0(name[part],'/MD/stabilisation/eqv_',name[part],'.conf'),row.names = F,col.names = F,quote = F)
  
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  #first MD
  df_tcl[1,1]<-paste0('structure protein/ionized_',name[part],'.psf',
                      '\ncoordinates pdb/eqv_',name[part],'.coor',
                      '\ntemperature 310',
                      '\n#############################################################',
                      '\n## SIMULATION PARAMETERS ##',
                      '\n#############################################################',
                      '\nparaTypeCharmm on',
                      '\nparameters ',part_start,'start/toppar/par_all36_prot.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_carb.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_cgenff.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_lipid.prm',
                      '\nparameters ',part_start,'start/toppar/par_all36_na.prm',
                      '\nparameters ',part_start,'start/toppar/toppar_water_ions_namd.str',
                      '\nexclude		scaled1-4',
                      '\n1-4scaling	1.0',
                      '\nswitching	on',
                      '\nswitchdist	8.0',
                      '\ncutoff		12.0',
                      '\npairlistdist	13.5',
                      '\nstepspercycle	20',
                      '\nPME		on',
                      '\nPMETolerance 	0.000001',
                      '\nPMEGridSizeX ',round(x_max-x_min,digits = 0),
                      '\nPMEGridSizeY ',round(y_max-y_min,digits = 0),	
                      '\nPMEGridSizeZ ',round(z_max-z_min,digits = 0),
                      '\nfullElectFrequency 4',
                      '\ncellBasisVector1 ',round(x_max-x_min,digits = 0),' 0.0 0.0',
                      '\ncellBasisVector2 0.0 ',round(y_max-y_min,digits = 0),' 0.0',
                      '\ncellBasisVector3 0.0 0.0 ',round(z_max-z_min,digits = 0),
                      '\ncellOrigin ',x_mean,' ',y_mean,' ',z_mean,
                      '\nwrapAll on',
                      '\ntimestep 1.0',
                      '\nnumsteps 10000000',
                      '\nseed 322223322', '\nuseGroupPressure yes', '\nuseFlexibleCell yes', '\nuseConstantRatio yes',
                      '\nLangevinPiston on', 
                      '\nLangevinPistonTarget 1.01325',
                      '\nLangevinPistonPeriod 200',
                      '\nLangevinPistonDecay 100',
                      '\nLangevinPistonTemp 310',
                      '\noutputenergies	1000',
                      '\ndcdfreq		1000',
                      '\ndcdfile quench/quench_',name[part],'_1.dcd',
                      '\nbinaryoutput no',
                      '\noutputname pdb/quench_',name[part],'_1')
  
  write.table(df_tcl,paste0(name[part],'/MD/stabilisation/quench_',name[part],'_1.conf'),row.names = F,col.names = F,quote = F)
  
  
  
  #MD 2-5
  
  for (j in 2:num_din) { 
    
    df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
    
    df_tcl[1,1]<-paste0('structure protein/ionized_',name[part],'.psf',
                        '\ncoordinates pdb/quench_',name[part],'_',(j-1),'.coor',
                        '\ntemperature 310',
                        '\n#############################################################',
                        '\n## SIMULATION PARAMETERS ##',
                        '\n#############################################################',
                        '\nparaTypeCharmm on',
                        '\nparameters ',part_start,'start/toppar/par_all36_prot.prm',
                        '\nparameters ',part_start,'start/toppar/par_all36_carb.prm',
                        '\nparameters ',part_start,'start/toppar/par_all36_cgenff.prm',
                        '\nparameters ',part_start,'start/toppar/par_all36_lipid.prm',
                        '\nparameters ',part_start,'start/toppar/par_all36_na.prm',
                        '\nparameters ',part_start,'start/toppar/toppar_water_ions_namd.str',
                        '\nexclude		scaled1-4',
                        '\n1-4scaling	1.0',
                        '\nswitching	on',
                        '\nswitchdist	8.0',
                        '\ncutoff		12.0',
                        '\npairlistdist	13.5',
                        '\nstepspercycle	20',
                        '\nPME		on',
                        '\nPMEGridSizeX ',round(x_max-x_min,digits = 0),
                        '\nPMEGridSizeY ',round(y_max-y_min,digits = 0),	
                        '\nPMEGridSizeZ ',round(z_max-z_min,digits = 0),
                        '\nfullElectFrequency 4',
                        '\ncellBasisVector1 ',round(x_max-x_min,digits = 0),' 0.0 0.0',
                        '\ncellBasisVector2 0.0 ',round(y_max-y_min,digits = 0),' 0.0',
                        '\ncellBasisVector3 0.0 0.0 ',round(z_max-z_min,digits = 0),
                        '\ncellOrigin ',x_mean,' ',y_mean,' ',z_mean,
                        '\nwrapAll on',
                        '\ntimestep 1.0',
                        '\nnumsteps 1000000',
                        '\nseed 322223322', '\nuseGroupPressure yes', '\nuseFlexibleCell yes', '\nuseConstantRatio yes',
                        '\nLangevinPiston on', 
                        '\nLangevinPistonTarget 1.01325',
                        '\nLangevinPistonPeriod 200',
                        '\nLangevinPistonDecay 100',
                        '\nLangevinPistonTemp 310',
                        '\noutputenergies	10000',
                        '\ndcdfreq		10000', 
                        '\ndcdfile quench/quench_',name[part],'_',j,'.dcd',
                        '\nbinaryoutput no',
                        '\noutputname pdb/quench_',name[part],'_',j)
    
    write.table(df_tcl,paste0(name[part],'/MD/stabilisation/quench_',name[part],'_',j,'.conf'),row.names = F,col.names = F,quote = F)
  }
}
#print MD script
df_conf<-data.frame(matrix(nrow = length(name),ncol = 1))
for (part in 1:length(name)) {
  df_conf[part,1]<-paste0(c(paste0("cd ",part_start,name[part],"/MD/stabilisation/\n"),
                            paste0(v_namd," " , c(paste0(" min_",name[part],".conf > min_",name[part],".out\n"),
                                     paste0(" heat_",name[part],".conf > heat_",name[part],".out\n"),
                                     paste0(" eqv_",name[part],".conf > eqv_",name[part],".out\n"),
                                     paste0(" quench_",name[part],"_",1:num_din,".conf > quench_",name[part],"_",1:num_din,".out\n")))),collapse = "")
}
write.table(df_conf,paste0(part_start,"r_scripts/namd_script.txt"),row.names = F,col.names = F,quote = F)
