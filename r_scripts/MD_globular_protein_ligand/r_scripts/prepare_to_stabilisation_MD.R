part_name = commandArgs(trailingOnly=TRUE)
# namd run script 
v_namd<-"namd run script"
#quantity of 25 ns MD simulations, change it to Modify length of MD simulation 

num_din<-1
library(dplyr)
library(bio3d)
setwd(part_name)

df_complex<-read.csv("start/complex.csv",stringsAsFactors = F)

df_complex<-unique(df_complex)
#df_replase<-read.csv("start/resname_mutate.csv",stringsAsFactors = F)
#df_complex<-left_join(df_complex,df_replase,by=c("ligand"="ligand_name"))
part<-1
for (part in 1:nrow(df_complex)) {
  print(paste0(df_complex$name.x[part]))
  if(!dir.exists(paste0("prepared_structures"))){dir.create(paste0("prepared_structures"))}
  if(!dir.exists(paste0("prepared_structures/start_structure"))){dir.create(paste0("prepared_structures/start_structure"))}
  if(!dir.exists(paste0("prepared_structures/prepared"))){dir.create(paste0("prepared_structures/prepared"))}
  if(!dir.exists(paste0("prepared_structures/chains"))){dir.create(paste0("prepared_structures/chains"))}
  if(!dir.exists(paste0("prepared_structures/tcl"))){dir.create(paste0("prepared_structures/tcl"))}
  if(!dir.exists(paste0("prepared_structures/complex_structure"))){dir.create(paste0("prepared_structures/complex_structure"))}
  if(!dir.exists(paste0("prepared_structures/fin_complex_structure"))){dir.create(paste0("prepared_structures/fin_complex_structure"))}
  if(!dir.exists(paste0("MD/"))){dir.create(paste0("MD/"))}
  if(!dir.exists(paste0("MD/",df_complex$number[part]))){dir.create(paste0("MD/",df_complex$number[part]))}
  if(!dir.exists(paste0("MD/",df_complex$number[part]))){dir.create(paste0("MD/",df_complex$number[part]))}
  #  system(command = paste0("cp -r ",part_name,"start/toppar/ ",part_name,"MD/",df_complex$number[part]),ignore.stdout=T,wait = T)
  #prepare psf and pdb parts of complex
  pdb_ligand<-read.pdb(paste0("start/ligands_center/",df_complex$name.x[part]))
  pdb_receptor<-read.pdb(paste0("start/receptor/",df_complex$receptor[part],".pdb"))
  pdb_ligand$atom$chain<-"Z"
  pdb_ligand$atom$elesy<-"Z"
  
  v_chain_ligand<-unique(pdb_ligand$atom$chain)
  v_chain_receptor<-unique(pdb_receptor$atom$chain)
  
  v_chain<-c(v_chain_receptor,v_chain_ligand)
  df_chain<-data.frame(matrix(nrow = length(v_chain),ncol = 2))
  colnames(df_chain)<-c("models","chain")
  df_chain$chain<-v_chain
  df_chain$models<-df_complex$name.x[part]
  
  df_chain_test<-df_chain%>%filter(is.na(chain))
  df_chain<-df_chain%>%filter(!is.na(chain))
  write.pdb(pdb_receptor,paste0("prepared_structures/start_structure/",df_chain$models[1],"_",df_chain$chain[1],".pdb"))
  write.pdb(pdb_ligand,paste0("prepared_structures/start_structure/",df_chain$models[2],"_",df_chain$chain[2],".pdb"))
  write.csv(df_chain,paste0("prepared_structures/chains/",df_complex$name.x[part],".csv"),row.names = F)
}
v_chains<-list.files(paste0("prepared_structures/chains/"))
df_chains<-read.csv(paste0("prepared_structures/chains/",v_chains[1]),stringsAsFactors = F)
i<-1
if(length(v_chains)>1){
  for (i in 2:length(v_chains)) {
    df_chains_add<-read.csv(paste0("prepared_structures/chains/",v_chains[i]),stringsAsFactors = F)
    df_chains<-rbind(df_chains,df_chains_add)
  }
}
df_chains<-left_join(df_complex,df_chains,by=c("name.x"="models"))
i<-1
for (i in 1:nrow(df_chains)) {
  df_psfgen<-data.frame(matrix(ncol = 1,nrow = 1))
  df_psfgen[1,1]<-paste0('cd ',part_name,'prepared_structures/\n',
                         'mol delete all\n',  
                         'package require psfgen\n',  
                         'resetpsf\n',  
                         'topology ',part_name,'start/toppar/toppar_water_ions.str\n',
                         'topology ',part_name,'start/toppar/toppar_ions_won.str\n',
                         'topology ',part_name,'start/toppar/toppar_dum_noble_gases.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_synthetic_polymer_patch.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_synthetic_polymer.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_prot_retinol.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_prot_na_combined.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_prot_modify_res.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_prot_model.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_prot_heme.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_prot_fluoro_alkanes.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_prot_c36m_d_aminoacids.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_prot_arg0.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_polymer_solvent.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_na_rna_modified.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_nano_lig_patch.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_nano_lig.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_na_nad_ppi.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_moreions.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_yeast.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_tag.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_sphingo.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_prot.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_mycobacterial.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_model.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_miscellaneous.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_lps.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_lnp.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_inositol.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_hmmm.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_ether.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_detergent.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_dag.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_cholesterol.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_cardiolipin.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_bacterial.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_lipid_archaeal.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_label_spin.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_label_fluorophore.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_carb_imlab.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_carb_glycopeptide.str\n',
                         'topology ',part_name,'start/toppar/toppar_all36_carb_glycolipid.str\n',
                         'topology ',part_name,'start/toppar/top_interface.rtf\n',
                         'topology ',part_name,'start/toppar/top_all36_prot.rtf\n',
                         'topology ',part_name,'start/toppar/top_all36_na.rtf\n',
                         'topology ',part_name,'start/toppar/top_all36_lipid.rtf\n',
                         'topology ',part_name,'start/toppar/top_all36_cgenff.rtf\n',
                         'topology ',part_name,'start/toppar/top_all36_carb.rtf\n',
                         'topology ',part_name,'start/toppar/tip216.crd\n',
                         'topology ',part_name,'start/toppar/cam.str\n')
  #  if (!is.na(df_chains$force_field[i])){  
  #  df_psfgen[2,1]<-paste0('topology ',part_name,'start/toppar/',df_chains$ligand[i],'\n')
  #  }
  df_psfgen[3,1]<-paste0('pdbalias residue HIS HSE\n',  
                         'pdbalias atom ILE CD1 CD\n')
  b<-paste0(strsplit(df_chains$ligand[i],split = "",fixed = T)[[1]][3:4],collapse = "")
  if ((df_chains$chain[i]=="Z")){  
    pdb<-read.pdb(paste0(part_name,'prepared_structures/start_structure/',df_chains$name.x[i],"_",df_chains$chain[i],'.pdb'))
    a<-unique(pdb$atom$resid)
    if(b=="CL"){
      df_psfgen[4,1]<-paste0('pdbalias residue ',a,' ', df_chains$ligand[i],2,'\n')
    }else{
      df_psfgen[4,1]<-paste0('pdbalias residue ',a,' ', df_chains$ligand[i],'\n')
    }
  }
  
  df_psfgen[5,1]<-paste0('segment ',df_chains$chain[i],' { pdb start_structure/',df_chains$name.x[i],"_",df_chains$chain[i],'.pdb\n',
                         '}',
                         '\ncoordpdb start_structure/',df_chains$name.x[i],"_",df_chains$chain[i],'.pdb ',df_chains$chain[i],'\n',
                         '\nregenerate angles dihedrals\n', 
                         '\nguesscoord\n',
                         
                         'writepdb prepared/',df_chains$name.x[i],"_",df_chains$chain[i],'.pdb\n',
                         'writepsf prepared/',df_chains$name.x[i],"_",df_chains$chain[i],'.psf\n',
                         'mol delete all \n \n \n exit now')
  write.table(df_psfgen,paste0('prepared_structures/tcl/psfgen_',df_chains$name.x[i],"_",df_chains$chain[i],'.tcl'),na = "\n",col.names = F,row.names = F,quote = F)
  system(command = paste0("vmd -dispdev text -e ",part_name,'prepared_structures/tcl/psfgen_',df_chains$name.x[i],"_",df_chains$chain[i],'.tcl'),ignore.stdout=T,wait = T) 
}
#df_complex<-read.csv("start/complex.csv",stringsAsFactors = F)
df_chains<-df_chains%>%mutate(structure_name=name.x)
v_complex<-unique(df_chains$structure_name)

i<-1

for (i in 1:length(v_complex)) {
  df_chain_TEMP<-df_chains%>%filter(structure_name==v_complex[i])
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


#Energy counting
i<-1
for (i in 1:length(v_complex)) {
  df_chain_TEMP<-df_chains%>%filter(structure_name==v_complex[i])
  df_chain_TEMP<-df_chain_TEMP%>%mutate(psf_merge=paste0(structure_name,"_",chain))
  df_chain_TEMP<-unique(df_chain_TEMP)
  
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  v_topology<-paste0(' -par ',part_name,'start/toppar/toppar_all36_synthetic_polymer_patch.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_synthetic_polymer.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_prot_retinol.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_prot_na_combined.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_prot_modify_res.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_prot_model.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_prot_heme.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_prot_fluoro_alkanes.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_prot_arg0.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_polymer_solvent.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_na_rna_modified.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_nano_lig_patch.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_nano_lig.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_na_nad_ppi.str ',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_tag.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_prot.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_mycobacterial.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_model.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_miscellaneous.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_lps.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_lnp.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_inositol.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_hmmm.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_ether.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_detergent.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_dag.str ',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_cardiolipin.str ',
                     ' -par ',part_name,'start/toppar/toppar_all36_lipid_archaeal.str ',
                     ' -par ',part_name,'start/toppar/toppar_all36_label_fluorophore.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_carb_imlab.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_carb_glycopeptide.str',
                     ' -par ',part_name,'start/toppar/toppar_all36_carb_glycolipid.str',
                     ' -par ',part_name,'start/toppar/par_interface.prm',
                     ' -par ',part_name,'start/toppar/par_all36_na.prm',
                     ' -par ',part_name,'start/toppar/par_all36m_prot.prm',
                     ' -par ',part_name,'start/toppar/par_all36_lipid.prm',
                     ' -par ',part_name,'start/toppar/par_all36_cgenff.prm',
                     ' -par ',part_name,'start/toppar/par_all36_carb.prm')
  df_tcl[1,1]<-paste0('cd ', part_name,'prepared_structures/\npackage require namdenergy')
  df_tcl[1,2]<-paste0('mol new {complex_structure/',v_complex[i],'.psf} type {psf}')
  df_tcl[1,3]<-paste0('mol addfile {complex_structure/',v_complex[i],'.pdb} type {pdb} waitfor all')
  df_tcl[1,4]<-paste0('set sel2 [atomselect top "protein"]')
  df_tcl[1,5]<-paste0('set sel1 [atomselect top "not (protein or water or ions)"]')
  df_tcl[1,6]<-paste0('\nnamdenergy -sel $sel1 $sel2 -vdw -elec -nonb -cutoff 12 -skip 0 -ofile Energy/interactions_',v_complex[i],'.txt -switch 10 -exe ',part_name,'programs/NAMD_2.14_Linux-x86_64-multicore/namd2 ',v_topology)
  
  #  df_tcl[1,7]<-paste0('\nnamdenergy -sel $sel2  -bond -angl -dihe -impr -conf -vdw -elec -nonb -all -cutoff 12 -skip 0 -ofile din/Energy/protein_',based_name[name],'.txt -switch 10 -exe /home/nastia/NAMD_2.14_Linux-x86_64-multicore/namd2 ',v_topology)
  df_tcl[1,7]<-'mol delete all'
  df_tcl[1,8]<-'\n\nexit now'
  write.table(df_tcl,file =paste0('prepared_structures/tcl/Energy_',v_complex[i],'.tcl'),sep = '\n',na = '' ,row.names = F,col.names = F,quote = F)
  print(paste0('vmd -dispdev text -e ',part_name,'prepared_structures/tcl/Energy_',v_complex[i],'.tcl'))
  system(command = paste0('vmd -dispdev text -e ',part_name,'prepared_structures/tcl/Energy_',v_complex[i],'.tcl'),ignore.stdout=T,wait = T) 
}
