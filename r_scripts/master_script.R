part_start<-"/home/nastia/projects/current/stabilisation_complex-receptor-ligand/"
setwd(part_start)

v_list_proteins<-list.files("start/sequence/")
a<-c()
for (i in 1:length(v_list_proteins)) {
  b<-strsplit(v_list_proteins[i],split = ".",fixed = T)[[1]][1]
  a<-c(a,b)
}
v_list_proteins<-a
for (i in 1:length(v_list_proteins)) {
  if(!dir.exists(v_list_proteins[i])){dir.create(v_list_proteins[i])}
}
#structure_prediction
for (i in 1:length(v_list_proteins)) {
  print(v_list_proteins[i])
  if (!dir.exists(paste0(part_start,v_list_proteins[i]))){dir.create(paste0(part_start,v_list_proteins[i]))}
  if (!dir.exists(paste0(part_start,v_list_proteins[i],"/structure_prediction/"))){dir.create(paste0(part_start,v_list_proteins[i],"/structure_prediction/"))}
  if (!dir.exists(paste0(part_start,v_list_proteins[i],"/structure_prediction/r_scripts/"))){dir.create(paste0(part_start,v_list_proteins[i],"/structure_prediction/r_scripts/"))}
  system(command = paste0("cp -r ",part_start,"start/structure_prediction/pdb/",v_list_proteins[i],"/ ",part_start,v_list_proteins[i],"/structure_prediction/pdb/"),ignore.stdout=T,wait = T)
  system(command = paste0("cp ",part_start,"start/structure_prediction/sequence/",v_list_proteins[i],".fasta ",part_start,v_list_proteins[i],"/structure_prediction/seqs.fasta"),ignore.stdout=T,wait = T)
  system(command = paste0("cp ",part_start,"start/sequence/",v_list_proteins[i],".fasta ",part_start,v_list_proteins[i],"/structure_prediction/protein_sequence.fasta"),ignore.stdout=T,wait = T)
  if(!dir.exists("predicted/")){dir.create("predicted/")}
  if(!dir.exists(paste0("predicted/",v_list_proteins[i]))){dir.create(paste0("predicted/",v_list_proteins[i]))}
  system(command = paste0("cp -r ",part_start,"r_scripts/structure_prediction/ ",part_start,v_list_proteins[i],"/structure_prediction/r_scripts/ "),ignore.stdout=T,wait = T)
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/structure_prediction.R ",part_start,v_list_proteins[i],"/"),ignore.stdout=T,wait = T)
  system(command = paste0("cp -r ",part_start,v_list_proteins[i],"/structure_prediction/pdb_fin/ ",part_start,"output/predicted/",v_list_proteins[i],"/"),ignore.stdout=T,wait = T)
  
}
#prepare files for MD based on chousen structures 
#predicted structures
#and sequense predict structure of protein using robetta server structure_prediction
i<-1
for (i in 1:length(v_list_proteins)) {
  print(v_list_proteins[i])
  if(!dir.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein/"))){dir.create(paste0(part_start,v_list_proteins[i],"/MD_globular_protein/"))}
  if(!dir.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein/programs/"))){dir.create(paste0(part_start,v_list_proteins[i],"/MD_globular_protein/programs/"))}
  system(command = paste0("cp -r ",part_start,"programs/NAMD_2.14_Linux-x86_64-multicore.tar.gz ",part_start,v_list_proteins[i],"/MD_globular_protein/programs/"),ignore.stdout=T,wait = T)
  
  system(command = paste0("cp -r ",part_start,"programs/la1.0 ",part_start,v_list_proteins[i],"/MD_globular_protein/programs/"),ignore.stdout=T,wait = T)
  system(command = paste0("cp -r ",part_start,"programs/orient ",part_start,v_list_proteins[i],"/MD_globular_protein/programs/"),ignore.stdout=T,wait = T)
  
#  system(command = paste0("tar -xvzf ",part_start,v_list_proteins[i],"/MD_globular_protein/programs/NAMD_2.14_Linux-x86_64-multicore.tar.gz " ,part_start,v_list_proteins[i],"/MD_globular_protein/programs/"),ignore.stdout=T,wait = T)
  
  system(command = paste0("rm ",part_start,"r_scripts/MD_globular_protein/start/structure/*"),ignore.stdout=T,wait = T)  
  system(command = paste0("cp ",part_start,"start/MD_stabilisation/",v_list_proteins[i],".pdb ",part_start,"r_scripts/MD_globular_protein/start/structure/",v_list_proteins[i],".pdb"),ignore.stdout=T,wait = T)
  system(command = paste0("cp -r ",part_start,"r_scripts/MD_globular_protein/ ",part_start,v_list_proteins[i],"/"),ignore.stdout=T,wait = T)
  system(command = paste0("cp -r ",part_start,"start/toppar/ ",part_start,v_list_proteins[i],"/MD_globular_protein/start/toppar/"),ignore.stdout=T,wait = T)
}
#MD stabilisation
i<-1
for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein/")
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein/r_scripts/prepare_to_stabilisation_MD.R ",part_name),ignore.stdout=T,wait = T)
  
#  print(paste0("run namd from ",part_name,"r_scripts/namd_script.txt"))
}
print(paste0("run namd from ",part_start,v_list_proteins,"/MD_globular_protein/r_scripts/namd_script.txt"))
i<-1
for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein/")
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein/r_scripts/prepare_tcl_din.R ",part_name),ignore.stdout=T,wait = T)
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein/r_scripts/second_stucture_compare.R ",part_name),ignore.stdout=T,wait = T)
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein/r_scripts/Ramachadran.R ",part_name),ignore.stdout=T,wait = T)
  #make plot
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein/r_scripts/make_plots_RMSD_RMSF.R ",part_name),ignore.stdout=T,wait = T)
}

for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein/")
  system(command = paste0("cp -r ",part_start,"r_scripts/MD_globular_protein/ ",part_start,v_list_proteins[i],"/"),ignore.stdout=T,wait = T)
  system(command = paste0("cp -r ",part_start,v_list_proteins[i],"/structure_prediction/pdb_fin/ ",part_start,"output/predicted/",v_list_proteins[i],"/"),ignore.stdout=T,wait = T)
}

#docking first 
#prepare ligands for first docking
v_ligands<-list.files("start/docking/docking_first/ligand_start/")
for (i in 1:length(v_ligands)) {
  a<-strsplit(v_ligands[i],split = ".",fixed = T)[[1]][1]
  system(command = paste0("obabel ",part_start,"start/docking/docking_first/ligand_start/",a, ".pdb -O ",part_start,"start/docking/docking_first/ligand/",a, ".pdbqt"),ignore.stdout=T,wait = T)
}
i<-1
for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,v_list_proteins[i],"/docking/")
  #copying sctipts for docking
  system(command = paste0("cp -r ",part_start,"r_scripts/docking/ ",part_start,v_list_proteins[i],"/"),ignore.stdout=T,wait = T)
  #create directories
  if(!dir.exists(paste0(part_name,"docking_first/"))){dir.create(paste0(part_name,"docking_first/"))}
  if(!dir.exists(paste0(part_name,"docking_first/receptor_start/"))){dir.create(paste0(part_name,"docking_first/receptor_start/"))}
  if(!dir.exists(paste0(part_name,"docking_first/receptor/"))){dir.create(paste0(part_name,"docking_first/receptor/"))}
  #copying ligand's structures for docking
  system(command = paste0("cp -r ",part_start,"start/docking/docking_first/ligand/ ",part_name,"/docking_first/"),ignore.stdout=T,wait = T)
  #copying receptor structure for docking
  system(command = paste0("cp -r ",part_start,"output/stabilisation/fin_structure/",v_list_proteins[i],".pdb ",part_name,"/docking_first/receptor_start/start.pdb"),ignore.stdout=T,wait = T)
  #copying amino acids for active center for docking
  system(command = paste0("cp ",part_start,"start/docking/active_center/",v_list_proteins[i],"/active_center.csv ",part_name,"/docking_first/active_center.csv"),ignore.stdout=T,wait = T)
  #copying programs for docking
  system(command = paste0("cp -r ",part_start,"programs/ ",part_name,"docking_first/"),ignore.stdout=T,wait = T)
}
i<-1
for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,v_list_proteins[i],"/docking/")
  system(command = paste0("Rscript --vanilla  ",part_name,"r_scripts/prepare_first_docking_main.R ",part_name),ignore.stdout=T,wait = T)
}
i<-1
for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,v_list_proteins[i],"/docking/")
  system(command = paste0("Rscript --vanilla  ",part_name,"r_scripts/first_docking_start_analysis.R ",part_name),ignore.stdout=T,wait = T)
}
i<-2
for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,v_list_proteins[i],"/docking/")
  part_scriprs<-paste0(part_start,"r_scripts/docking/r_scripts/")
  part_analysis<-paste0(part_name,"docking_first/")
  system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_pre_analysis.R ",part_analysis),ignore.stdout=T,wait = T)
  system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_group_structure.R ",part_analysis),ignore.stdout=T,wait = T)
  system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_interactions.R ",part_analysis),ignore.stdout=T,wait = T)
  system(command = paste0("Rscript --vanilla  ",part_scriprs,"analysis.R ",part_analysis,""),ignore.stdout=T,wait = T)
}


#MD receptor-ligand stabilisation
i<-1
j<-1
for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/")
  print(v_list_proteins[i])
  system(command = paste0("cp -r ",part_start,"programs/NAMD_2.14_Linux-x86_64-multicore.tar.gz ",part_start,v_list_proteins[i],"/MD_globular_protein/programs/"),ignore.stdout=T,wait = T)
  if(!dir.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/"))){dir.create(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/"))}
  if(!dir.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/programs/"))){dir.create(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/programs/"))}
  if(!dir.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/"))){dir.create(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/"))}
  if(!dir.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/ligands/"))){dir.create(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/ligands/"))}
  if(!dir.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/receptor/"))){dir.create(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/receptor/"))}
 
  system(command = paste0("cp -r ",part_start,v_list_proteins[i],"/docking/docking_first/receptor_start/start.pdb ", part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/receptor/start.pdb"))
  system(command = paste0("cp -r ",part_start,"programs/NAMD_2.14_Linux-x86_64-multicore.tar.gz ",part_start,v_list_proteins[i],"/MD_globular_protein_ligand/programs/"),ignore.stdout=T,wait = T)
#  system(command = paste0("tar -xvzf ",part_start,v_list_proteins[i],"/MD_globular_protein_ligand/programs/NAMD_2.14_Linux-x86_64-multicore.tar.gz ",part_start,v_list_proteins[i],"/MD_globular_protein_ligand/programs/"),ignore.stdout=T,wait = T)


  v_ligands<-list.files(paste0(part_start,v_list_proteins[i],"/docking/docking_first/din/str_fin/"))
  if (length(v_ligands)>0){
    for (j in 1:length(v_ligands)) {
      system(command = paste0("cp -r ",part_start,v_list_proteins[i],"/docking/docking_first/din/str_fin/",v_ligands[j]," ",part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/ligands/"),ignore.stdout=T,wait = T)
    }
  }
  if(file.exists(paste0(part_start,v_list_proteins[i],"/docking/docking_first/receptor_start/start.pdb"))){
    system(command = paste0("cp -r ",part_start,v_list_proteins[i],"/docking/docking_first/receptor_start/start.pdb ", part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/receptor/start.pdb"))
  }

  system(command = paste0("cp -r ",part_start,"r_scripts/MD_globular_protein_ligand/ ",part_start,v_list_proteins[i],"/"),ignore.stdout=T,wait = T)
  system(command = paste0("cp -r ",part_start,"start/toppar/ ",part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/toppar/"),ignore.stdout=T,wait = T)
  system(command = paste0("cp -r ",part_start,"start/receptor_ligand_MD/resname_mutate.csv ",part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/"),ignore.stdout=T,wait = T)
  system(command = paste0("cp -r ",part_start,"r_scripts/MD_globular_protein_ligand/ ",part_start,v_list_proteins[i],"/"),ignore.stdout=T,wait = T)
  
}
i<-1

for (i in 1:length(v_list_proteins)) {
  if(file.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/receptor/start.pdb"))){
    part_prepare<-paste0(part_start,v_list_proteins[i])
    part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/")
    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein_ligand/r_scripts/convert_docking_MD_receptor_ligand.R ",part_prepare),ignore.stdout=T,wait = T)
    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein_ligand/r_scripts/prepare_to_stabilisation_MD.R ",part_name),ignore.stdout=T,wait = T)
    system(command = paste0("vmd -dispdev text -e ",part_start,"r_scripts/MD_globular_protein_ligand/r_scripts/prepare_MD.tcl "),ignore.stdout=T,wait = T) 
    #namd
    print(paste0("#run namd2 from "))
    print(paste0(part_start,"r_scripts/namd_script.txt"))
  }
}
i<-1
j<-1
for (i in 1:length(v_list_proteins)) {
  part_prepare<-paste0(part_start,"r_scripts/MD_globular_protein_ligand/")
  part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/") 
  if (file.exists(paste0(part_name,"complex.csv"))){
    df_complex<-read.csv(paste0(part_name,"complex.csv"),stringsAsFactors = F)
    for (j in 1:nrow(df_complex)) {
      if (length(list.files(paste0(part_name,df_complex$complex_name[j],"/MD/stabilisation/quench/")))>0){
        print(paste0(part_name,df_complex$complex_name[j]))
        part_analysis<-paste0(part_name,df_complex$complex_name[j],"/")    
        system(command = paste0("Rscript --vanilla  ",part_prepare,"r_scripts/prepare_tcl_din.R ",part_analysis),ignore.stdout=T,wait = T)
        system(command = paste0("Rscript --vanilla  ",part_prepare,"r_scripts/second_stucture_compare.R ",part_analysis),ignore.stdout=T,wait = T)
        system(command = paste0("Rscript --vanilla  ",part_prepare,"r_scripts/Ramachadran.R ",part_analysis),ignore.stdout=T,wait = T)
        #make plot
        system(command = paste0("Rscript --vanilla  ",part_prepare,"r_scripts/make_plots_RMSD_RMSF.R ",part_analysis),ignore.stdout=T,wait = T)
      }
    }
  }
}

i<-1
j<-1
for (i in 1:length(v_list_proteins)) {
  part_prepare<-paste0(part_start,"r_scripts/MD_globular_protein_ligand/")
  part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/") 
#  part_start<-part_name
  if (file.exists(paste0(part_name,"complex.csv"))){
    df_complex<-read.csv(paste0(part_name,"complex.csv"),stringsAsFactors = F)
    for (j in 1:nrow(df_complex)) {
      if (length(list.files(paste0(part_name,df_complex$complex_name[j],"/MD/stabilisation/quench/")))>0){
        print(paste0(part_name,df_complex$complex_name[j]))
        part_analysis<-paste0(part_name,df_complex$complex_name[j],"/")    
        system(command = paste0("Rscript --vanilla  ",part_prepare,"r_scripts/complex_structure_analysis.R ",part_analysis),ignore.stdout=T,wait = T)
      }
    }
  }
}
i<-1
for (i in 1:length(v_list_proteins)) {
  part_prepare<-paste0(part_start,"r_scripts/MD_globular_protein_ligand/")
  part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/") 
  system(command = paste0("Rscript --vanilla  ",part_prepare,"r_scripts/fin_structure_prepare.R ",part_name),ignore.stdout=T,wait = T)
  system(command = paste0("Rscript --vanilla  ",part_prepare,"r_scripts/atom_interactions.R ",part_name),ignore.stdout=T,wait = T)
  system(command = paste0("Rscript --vanilla  ",part_prepare,"r_scripts/make_structure_picture.R ",part_name),ignore.stdout=T,wait = T)
}
