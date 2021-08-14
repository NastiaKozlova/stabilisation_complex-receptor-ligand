part_name = commandArgs(trailingOnly=TRUE)
# namd run script 
v_namd<-"namd run script"
#quantity of 25 ns MD simulations, change it to Modify length of MD simulation 

num_din<-5
library(dplyr)
library(bio3d)
setwd(part_name)
df_complex<-read.csv("start/complex.csv",stringsAsFactors = F)
#df_complex<-df_complex%>%select(structure_name,resname_change)
df_complex<-unique(df_complex)
df_replase<-read.csv("start/resname_mutate.csv",stringsAsFactors = F)
df_complex<-left_join(df_complex,df_replase,by=c("ligand"="ligand_name"))
part<-1
for (part in 1:nrow(df_complex)) {
  print(paste0(df_complex$structure_name[part]))
  if(!dir.exists(paste0("prepared_structures"))){dir.create(paste0("prepared_structures"))}
  if(!dir.exists(paste0("prepared_structures/start_structure"))){dir.create(paste0("prepared_structures/start_structure"))}
  if(!dir.exists(paste0("prepared_structures/prepared"))){dir.create(paste0("prepared_structures/prepared"))}
  if(!dir.exists(paste0("prepared_structures/chains"))){dir.create(paste0("prepared_structures/chains"))}
  if(!dir.exists(paste0("prepared_structures/tcl"))){dir.create(paste0("prepared_structures/tcl"))}
  if(!dir.exists(paste0("prepared_structures/complex_structure"))){dir.create(paste0("prepared_structures/complex_structure"))}
  if(!dir.exists(paste0("prepared_structures/fin_complex_structure"))){dir.create(paste0("prepared_structures/fin_complex_structure"))}
  #prepare psf and pdb parts of complex
  pdb_ligand<-read.pdb(paste0("start/ligands/",df_complex$structure_name[part]))
  pdb_receptor<-read.pdb(paste0("start/receptor/",df_complex$receptor[part],".pdb"))
  pdb_ligand$atom$chain<-"Z"
  pdb_ligand$atom$elesy<-"Z"
#  pdb_ligand$atom$elesy<-""
  v_chain_ligand<-unique(pdb_ligand$atom$chain)
  v_chain_receptor<-unique(pdb_receptor$atom$chain)
  
  v_chain<-c(v_chain_receptor,v_chain_ligand)
  df_chain<-data.frame(matrix(nrow = length(v_chain),ncol = 2))
  colnames(df_chain)<-c("models","chain")
  df_chain$chain<-v_chain
  df_chain$models<-df_complex$structure_name[part]
  
  df_chain_test<-df_chain%>%filter(is.na(chain))
  df_chain<-df_chain%>%filter(!is.na(chain))
  write.pdb(pdb_receptor,paste0("prepared_structures/start_structure/",df_chain$models[1],"_",df_chain$chain[1],".pdb"))
  write.pdb(pdb_ligand,paste0("prepared_structures/start_structure/",df_chain$models[2],"_",df_chain$chain[2],".pdb"))
  write.csv(df_chain,paste0("prepared_structures/chains/",df_complex$structure_name[part],".csv"),row.names = F)
}
v_chains<-list.files(paste0("prepared_structures/chains/"))
df_chains<-read.csv(paste0("prepared_structures/chains/",v_chains[1]),stringsAsFactors = F)
#i<-2
if(length(v_chains)>1){
  for (i in 2:length(v_chains)) {
    df_chains_add<-read.csv(paste0("prepared_structures/chains/",v_chains[i]),stringsAsFactors = F)
    df_chains<-rbind(df_chains,df_chains_add)
  }
}
df_chains<-left_join(df_complex,df_chains,by=c("structure_name"="models"))
i<-1

for (i in 1:nrow(df_chains)) {
  df_psfgen<-data.frame(matrix(ncol = 1,nrow = 3))
  df_psfgen[1,1]<-paste0('cd ',part_name,'prepared_structures/\n',
                         'mol delete all\n',  
                         'package require psfgen\n',  
                         'resetpsf\n',  
                         'topology ',part_name,'start/toppar/top_all36_prot.rtf\n',  
                         'topology ',part_name,'start/toppar/toppar_water_ions_namd.str\n',  
                         'topology ',part_name,'start/toppar/top_all36_carb.rtf\n',  
                         'topology ',part_name,'start/toppar/top_all36_cgenff.rtf\n',  
                         'topology ',part_name,'start/toppar/top_all36_lipid.rtf\n',  
                         'topology ',part_name,'start/toppar/top_all36_na.rtf\n',  
                         'topology ',part_name,'start/toppar/top_all36_prot.rtf\n')
  if (!is.na(df_chains$force_field[i])){  
    df_psfgen[2,1]<-paste0('topology ',part_name,'start/toppar/',df_chains$force_field[i],'\n')
  }
  df_psfgen[3,1]<-paste0('pdbalias residue HIS HSE\n',  
                         'pdbalias atom ILE CD1 CD\n')
  if (!is.na(df_chains$resname_change[i])){  
    df_psfgen[4,1]<-paste0('pdbalias residue ',df_chains$pdb_name[i],' ', df_chains$charmm_name[i],'\n')}
  
  df_psfgen[5,1]<-paste0('segment ',df_chains$chain[i],' { pdb start_structure/',df_chains$structure_name[i],"_",df_chains$chain[i],'.pdb\n',
                         '}',
                         '\ncoordpdb start_structure/',df_chains$structure_name[i],"_",df_chains$chain[i],'.pdb ',df_chains$chain[i],'\n',  
                         'regenerate angles dihedrals\n', 
                         'guesscoord\n',
                         'writepdb prepared/',df_chains$structure_name[i],"_",df_chains$chain[i],'.pdb\n',
                         'writepsf prepared/',df_chains$structure_name[i],"_",df_chains$chain[i],'.psf\n',
                         'mol delete all \n \n \n exit now')
  write.table(df_psfgen,paste0('prepared_structures/tcl/psfgen_',df_chains$structure_name[i],"_",df_chains$chain[i],'.tcl'),na = "\n",col.names = F,row.names = F,quote = F)
  system(command = paste0("vmd -dispdev text -e ",part_name,'prepared_structures/tcl/psfgen_',df_chains$structure_name[i],"_",df_chains$chain[i],'.tcl'),ignore.stdout=T,wait = T) 
}
df_complex<-read.csv("start/complex.csv",stringsAsFactors = F)
v_complex<-unique(df_chains$complex)
i<-2
for (i in 1:length(v_complex)) {
  df_chain_TEMP<-df_chains%>%filter(complex_name==v_complex[i])
  df_chain_TEMP<-df_chain_TEMP%>%mutate(psf_merge=paste0(structure_name,"_",chain))
  df_chain_TEMP<-unique(df_chain_TEMP)

  df_tcl<-  c(paste0("cd ",part_name,"prepared_structures/\n",
                     "mol delete all\n",
                     "package require psfgen\n",
                     "resetpsf\n"),
              
              
              unique(paste0("readpsf prepared/",df_chain_TEMP$psf_merge,".psf \ncoordpdb prepared/",df_chain_TEMP$psf_merge,".pdb\n")),
              
              paste0("regenerate angles dihedrals\n",
                     "guesscoord\n",
                     "writepdb complex_structure/",v_complex[i],".pdb\n",
                     "writepsf complex_structure/",v_complex[i],".psf\n\n",
                     'mol delete all \n \n \n exit now')
  )
  df_tcl<-as.data.frame(df_tcl)
  write.table(df_tcl,paste0('prepared_structures/tcl/combine_',v_complex[i],'.tcl'),na = "\n",col.names = F,row.names = F,quote = F)
  system(command = paste0("vmd -dispdev text -e ",part_name,'prepared_structures/tcl/combine_',v_complex[i],'.tcl'),ignore.stdout=T,wait = T) 
}
part<-2
for (part in 1:length(v_complex)) {
  print(paste0(v_complex[part]))
  if(file.exists(paste0("prepared_structures/complex_structure/",v_complex[part],'.pdb'))){
    pdb_start<-read.pdb(paste0("prepared_structures/complex_structure/",v_complex[part],'.pdb'))
    df_pdb<-pdb_start$atom
    
    x_min<-min(df_pdb$x)-10
    y_min<-min(df_pdb$y)-10
    z_min<-min(df_pdb$z)-10
    
    x_max<-max(df_pdb$x)+10
    y_max<-max(df_pdb$y)+10
    z_max<-max(df_pdb$z)+10
    #prepare script to solvate and ionize protein pstructure
    df_psfgen<-data.frame(matrix(ncol = 1,nrow = 1))
    df_psfgen[1,1]<-paste0("cd ",part_name,"prepared_structures/\n",
                           'package require solvate \n','package require autoionize \n',
                           'solvate complex_structure/',v_complex[part],'.psf ',
                           'complex_structure/',v_complex[part],'.pdb',
                           ' -o complex_structure/solvate_',v_complex[part],' -b 1.5 -minmax {{',x_min,' ',y_min,' ',z_min,'} {',x_max,' ',y_max,' ',z_max,'}}\n',
                           'autoionize -psf complex_structure/solvate_',v_complex[part],'.psf -pdb complex_structure/solvate_',v_complex[part],'.pdb -sc 0.15 -o fin_complex_structure/ionized_',v_complex[part],
                           '\nmol delete all\n\nexit now')
    write.table(df_psfgen,paste0('prepared_structures/tcl/solvate_',v_complex[part],'.tcl'),col.names = F,row.names = F,quote = F)
    system(command = paste0("vmd -dispdev text -e ",part_name,'prepared_structures/tcl/solvate_',v_complex[part],'.tcl'),ignore.stdout=T,wait = T) 
  }
}
i<-2
for (i in 1:length(v_complex)) {
  if(!dir.exists(paste0(v_complex[i]))){dir.create(paste0(v_complex[i]))}
  if(!dir.exists(paste0(v_complex[i],"/MD"))){dir.create(paste0(v_complex[i],"/MD"))}
  if(!dir.exists(paste0(v_complex[i],"/MD/stabilisation/"))){dir.create(paste0(v_complex[i],"/MD/stabilisation/"))}
  if(!dir.exists(paste0(v_complex[i],"/MD/stabilisation/quench/"))){dir.create(paste0(v_complex[i],"/MD/stabilisation/quench/"))}
  if(!dir.exists(paste0(v_complex[i],"/MD/stabilisation/protein/"))){dir.create(paste0(v_complex[i],"/MD/stabilisation/protein/"))}
  if(!dir.exists(paste0(v_complex[i],"/MD/stabilisation/dcd/"))){dir.create(paste0(v_complex[i],"/MD/stabilisation/dcd/"))}
  if(!dir.exists(paste0(v_complex[i],"/MD/stabilisation/din/"))){dir.create(paste0(v_complex[i],"/MD/stabilisation/din/"))}
  if(!dir.exists(paste0(v_complex[i],"/MD/stabilisation/pdb/"))){dir.create(paste0(v_complex[i],"/MD/stabilisation/pdb/"))}
  
  if(file.exists(paste0(part_name,'prepared_structures/fin_complex_structure/ionized_',v_complex[i],".psf"))){
    system(command = paste0("cp ",part_name,'prepared_structures/fin_complex_structure/ionized_',v_complex[i],".psf ",part_name,v_complex[i],"/MD/stabilisation/protein/"),ignore.stdout=T,wait = T) 
    system(command = paste0("cp ",part_name,'prepared_structures/fin_complex_structure/ionized_',v_complex[i],".pdb ",part_name,v_complex[i],"/MD/stabilisation/protein/"),ignore.stdout=T,wait = T) 
    pdb_ionized<-read.pdb(paste0(v_complex[i],'/MD/stabilisation/protein/ionized_',v_complex[i],'.pdb'))
    df_pdb<-pdb_ionized$atom
    x_min<-min(df_pdb$x)
    y_min<-min(df_pdb$y)
    z_min<-min(df_pdb$z)
    
    x_max<-max(df_pdb$x)
    y_max<-max(df_pdb$y)
    z_max<-max(df_pdb$z)
    
    x_mean<-mean(df_pdb$x)
    y_mean<-mean(df_pdb$y)
    z_mean<-mean(df_pdb$z)
    #energy minimisation 
    df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
    df_tcl[1,1]<-paste0('structure protein/ionized_',v_complex[i],'.psf',
                        '\ncoordinates protein/ionized_',v_complex[i],'.pdb',
                        '\nset temperature 300',
                        '\nfirsttimestep 0',
                        '\n#############################################################',
                        '\n## SIMULATION PARAMETERS ##',
                        '\n#############################################################',
                        '\nparaTypeCharmm on',
                        '\nparameters ',part_name,'start/toppar/par_all36_carb.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_cgenff.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_lipid.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36m_prot.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_na.prm',
                        '\nparameters ',part_name,'start/toppar/toppar_water_ions_namd.str',
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
                        '\ndcdfile dcd/min_',v_complex[i],'.dcd',
                        '\nbinaryoutput no',
                        '\noutputname pdb/min_',v_complex[i])
    write.table(df_tcl,paste0(v_complex[i],'/MD/stabilisation/min_',v_complex[i],'.conf'),row.names = F,col.names = F,quote = F)
    
    df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
    #structure heating 
    df_tcl[1,1]<-paste0('structure protein/ionized_',v_complex[i],'.psf',
                        '\ncoordinates pdb/min_',v_complex[i],'.coor',
                        '\n#############################################################',
                        '\n## SIMULATION PARAMETERS ##',
                        '\n#############################################################',
                        '\nparaTypeCharmm on',
                        '\nparameters ',part_name,'start/toppar/par_all36_prot.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_carb.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_cgenff.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_lipid.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_na.prm',
                        '\nparameters ',part_name,'start/toppar/toppar_water_ions_namd.str',
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
                        '\ndcdfile dcd/head_',v_complex[i],'.dcd',
                        '\nbinaryoutput no',
                        '\noutputname pdb/head_',v_complex[i])
    write.table(df_tcl,paste0(v_complex[i],'/MD/stabilisation/heat_',v_complex[i],'.conf'),row.names = F,col.names = F,quote = F)  
    df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
    #first eqvilibration
    df_tcl[1,1]<-paste0('structure protein/ionized_',v_complex[i],'.psf',
                        '\ncoordinates pdb/head_',v_complex[i],'.coor',
                        '\n#############################################################',
                        '\n## SIMULATION PARAMETERS ##',
                        '\n#############################################################',
                        '\nparaTypeCharmm on',
                        '\nparameters ',part_name,'start/toppar/par_all36_prot.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_carb.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_cgenff.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_lipid.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_na.prm',
                        '\nparameters ',part_name,'start/toppar/toppar_water_ions_namd.str',
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
                        '\ndcdfile dcd/eqv_',v_complex[i],'.dcd',
                        '\nbinaryoutput no',
                        '\noutputname pdb/eqv_',v_complex[i])
    
    write.table(df_tcl,paste0(v_complex[i],'/MD/stabilisation/eqv_',v_complex[i],'.conf'),row.names = F,col.names = F,quote = F)
    
    df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
    #first MD
    df_tcl[1,1]<-paste0('structure protein/ionized_',v_complex[i],'.psf',
                        '\ncoordinates pdb/eqv_',v_complex[i],'.coor',
                        '\ntemperature 310',
                        '\n#############################################################',
                        '\n## SIMULATION PARAMETERS ##',
                        '\n#############################################################',
                        '\nparaTypeCharmm on',
                        '\nparameters ',part_name,'start/toppar/par_all36_prot.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_carb.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_cgenff.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_lipid.prm',
                        '\nparameters ',part_name,'start/toppar/par_all36_na.prm',
                        '\nparameters ',part_name,'start/toppar/toppar_water_ions_namd.str',
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
                        '\nnumsteps 25000000',
                        '\nseed 322223322', '\nuseGroupPressure yes', '\nuseFlexibleCell yes', '\nuseConstantRatio yes',
                        '\nLangevinPiston on', 
                        '\nLangevinPistonTarget 1.01325',
                        '\nLangevinPistonPeriod 200',
                        '\nLangevinPistonDecay 100',
                        '\nLangevinPistonTemp 310',
                        '\noutputenergies	10000',
                        '\ndcdfreq		10000',
                        '\ndcdfile quench/quench_',v_complex[i],'_1.dcd',
                        '\nbinaryoutput no',
                        '\noutputname pdb/quench_',v_complex[i],'_1')
    
    write.table(df_tcl,paste0(v_complex[i],'/MD/stabilisation/quench_',v_complex[i],'_1.conf'),row.names = F,col.names = F,quote = F)
    
    
    
    #MD 2-5
    
    for (j in 2:num_din) { 
      
      df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
      
      df_tcl[1,1]<-paste0('structure protein/ionized_',v_complex[i],'.psf',
                          '\ncoordinates pdb/quench_',v_complex[i],'_',(j-1),'.coor',
                          '\ntemperature 310',
                          '\n#############################################################',
                          '\n## SIMULATION PARAMETERS ##',
                          '\n#############################################################',
                          '\nparaTypeCharmm on',
                          '\nparameters ',part_name,'start/toppar/par_all36_prot.prm',
                          '\nparameters ',part_name,'start/toppar/par_all36_carb.prm',
                          '\nparameters ',part_name,'start/toppar/par_all36_cgenff.prm',
                          '\nparameters ',part_name,'start/toppar/par_all36_lipid.prm',
                          '\nparameters ',part_name,'start/toppar/par_all36_na.prm',
                          '\nparameters ',part_name,'start/toppar/toppar_water_ions_namd.str',
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
                          '\nnumsteps 25000000',
                          '\nseed 322223322', '\nuseGroupPressure yes', '\nuseFlexibleCell yes', '\nuseConstantRatio yes',
                          '\nLangevinPiston on', 
                          '\nLangevinPistonTarget 1.01325',
                          '\nLangevinPistonPeriod 200',
                          '\nLangevinPistonDecay 100',
                          '\nLangevinPistonTemp 310',
                          '\noutputenergies	10000',
                          '\ndcdfreq		10000', 
                          '\ndcdfile quench/quench_',v_complex[i],'_',j,'.dcd',
                          '\nbinaryoutput no',
                          '\noutputname pdb/quench_',v_complex[i],'_',j)
      
      write.table(df_tcl,paste0(v_complex[i],'/MD/stabilisation/quench_',v_complex[i],'_',j,'.conf'),row.names = F,col.names = F,quote = F)
    }
  }
}
#print MD script
df_conf<-data.frame(matrix(nrow = length(v_complex),ncol = 1))
for (part in 1:length(v_complex)) {
  if(file.exists(paste0(v_complex[part],'/MD/stabilisation/protein/ionized_',v_complex[part],'.pdb'))){
    df_conf[part,1]<-paste0(c(paste0("cd ",part_name,v_complex[part],"/MD/stabilisation/\n"),
                              paste0(v_namd," " , c(paste0(" min_",v_complex[part],".conf > min_",v_complex[part],".out\n"),
                                                    paste0(" heat_",v_complex[part],".conf > heat_",v_complex[part],".out\n"),
                                                    paste0(" eqv_",v_complex[part],".conf > eqv_",v_complex[part],".out\n"),
                                                    paste0(" quench_",v_complex[part],"_",1:num_din,".conf > quench_",v_complex[part],"_",1:num_din,".out\n")))),collapse = "")
    
  }
}
write.table(df_conf,paste0(part_name,"r_scripts/namd_script.txt"),row.names = F,col.names = F,quote = F,na = "")
