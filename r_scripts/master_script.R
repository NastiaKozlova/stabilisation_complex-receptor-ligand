part_start<-"path to stabilisation_complex-receptor-ligand/"
part_start<-paste0(getwd(),"/")
setwd(part_start)
#if you want don't count cout interactions of protein with protein surface surphase_conut<-F
surphase_conut<-T

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
i<-1
for (i in 1:length(v_list_proteins)) {
  print(v_list_proteins[i])
  part_name<-paste0(part_start,",",v_list_proteins[i])
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/structure_prediction.R ",part_name),ignore.stdout=T,wait = T)
  
}
#prepare files for MD based on chousen structures 
#predicted structures
#and sequense predict structure of protein using robetta server structure_prediction

for (i in 1:length(v_list_proteins)) {
  print(v_list_proteins[i])
  part_name<-paste0(part_start,",",v_list_proteins[i])
  print(part_name)
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein/r_scripts/copy_pre_MD_simulation_files.R ",part_name),ignore.stdout=T,wait = T)
}
#MD stabilisation
#prepare files to run MD simulation

for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein/")
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein/r_scripts/prepare_to_stabilisation_MD.R ",part_name),ignore.stdout=T,wait = T)
#  print(paste0("run namd from ",part_name,"r_scripts/namd_script.txt"))
}
#copy files to MD simulation
for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein/")
  system(command = paste0("cp -r ",part_start,"r_scripts/MD_globular_protein/ ",part_start,v_list_proteins[i],"/"),ignore.stdout=T,wait = T)
  system(command = paste0("cp -r ",part_start,v_list_proteins[i],"/structure_prediction/pdb_fin/ ",part_start,"output/predicted/",v_list_proteins[i],"/"),ignore.stdout=T,wait = T)
}
#generated scipts to run MD simulation are in this folder
#to continue experiment you should run MD simulation
print(paste0("run namd from ",part_start,v_list_proteins,"/MD_globular_protein/r_scripts/namd_script.txt"))

for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein/")
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein/r_scripts/prepare_tcl_din.R ",part_name),ignore.stdout=T,wait = T)
  part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein/")
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein/r_scripts/second_stucture_compare.R ",part_name),ignore.stdout=T,wait = T)
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein/r_scripts/Ramachadran.R ",part_name),ignore.stdout=T,wait = T)
    #make plot
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein/r_scripts/make_plots_RMSD_RMSF.R ",part_name),ignore.stdout=T,wait = T)
}
#copy final receptor structure to output folder
for (i in 1:length(v_list_proteins)) {
  if(!dir.exists(paste0(part_start,"output"))){dir.create(paste0(part_start,"output"))}
  if(!dir.exists(paste0(part_start,"output/stabilisation"))){dir.create(paste0(part_start,"output/stabilisation"))}
  system(command = paste0("cp ",part_start,v_list_proteins[i],"/MD_globular_protein/",v_list_proteins[i],"/MD/stabilisation/din/fin_structure/",v_list_proteins[i],".pdb",
                          " ",part_start,"output/stabilisation"),ignore.stdout=T,wait = T)
}

#docking first 
#prepare ligands for first docking

v_ligands<-list.files("start/docking/docking_first/ligand_start/")
if(!dir.exists(paste0(part_start,"start/docking/docking_first/ligand/"))){dir.create(paste0(part_start,"start/docking/docking_first/ligand/"))}
for (i in 1:length(v_ligands)) {
  a<-strsplit(v_ligands[i],split = ".",fixed = T)[[1]][1]
  system(command = paste0("obabel ",part_start,"start/docking/docking_first/ligand_start/",a, ".pdb -O ",part_start,"start/docking/docking_first/ligand/",a, ".pdbqt"),ignore.stdout=T,wait = T)
}

for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,",",v_list_proteins[i])
  #copying sctipts for docking
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/copy_files_for_docking.R ",part_name),ignore.stdout=T,wait = T)
}
i<-1
if(surphase_conut){
  #docking_main_surphase
  for (i in 1:length(v_list_proteins)) {
    part_protein<-paste0(part_start,",",v_list_proteins[i])
    #copying sctipts for docking
    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_main_surphase.R ",part_protein),ignore.stdout=T,wait = T)
  }
}else{
  for (i in 1:length(v_list_proteins)) {
    #prepare castade start_file
    part_protein<-paste0(part_start,",",v_list_proteins[i])
    #copying sctipts for docking
    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_main_center.R ",part_protein),ignore.stdout=T,wait = T)
    #sctipts for cascade docking
    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/cascade_docking_center.R ",part_protein),ignore.stdout=T,wait = T)
  }
}


#MD receptor-ligand stabilisation
#copy_files_for_MD_receptor_ligand
i<-1
for (i in 1:length(v_list_proteins)) {
  part_protein<-paste0(part_start,",",v_list_proteins[i])
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/copy_files_for_MD_receptor_ligand.R ",part_protein),ignore.stdout=T,wait = T)
  
#  copy_files_for_MD_receptor_ligand
#  part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/")
#  print(v_list_proteins[i])
#  system(command = paste0("cp -r ",part_start,"programs/NAMD_2.14_Linux-x86_64-multicore.tar.gz ",part_start,v_list_proteins[i],"/MD_globular_protein/programs/"),ignore.stdout=T,wait = T)
#  if(!dir.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/"))){dir.create(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/"))}
#  if(!dir.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/programs/"))){dir.create(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/programs/"))}
#  if(!dir.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/"))){dir.create(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/"))}
#  if(!dir.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/ligands/"))){dir.create(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/ligands/"))}
#  if(!dir.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/receptor/"))){dir.create(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/receptor/"))}
 
#  system(command = paste0("cp -r ",part_start,v_list_proteins[i],"/docking/docking_first/receptor_start/start.pdb ", part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/receptor/start.pdb"))
#  system(command = paste0("cp -r ",part_start,"programs/NAMD_2.14_Linux-x86_64-multicore.tar.gz ",part_start,v_list_proteins[i],"/MD_globular_protein_ligand/programs/"),ignore.stdout=T,wait = T)
#  system(command = paste0("tar -xvzf ",part_start,v_list_proteins[i],"/MD_globular_protein_ligand/programs/NAMD_2.14_Linux-x86_64-multicore.tar.gz ",part_start,v_list_proteins[i],"/MD_globular_protein_ligand/programs/"),ignore.stdout=T,wait = T)


#  v_ligands<-list.files(paste0(part_start,v_list_proteins[i],"/docking/docking_first/din/structure_merged/"))
#  if (length(v_ligands)>0){
#    for (j in 1:length(v_ligands)) {
#      system(command = paste0("cp -r ",part_start,v_list_proteins[i],"/docking/docking_first/din/structure_merged/",v_ligands[j]," ",part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/ligands/"),ignore.stdout=T,wait = T)
#    }
#  }
##  system(command = paste0("cp -r ",part_start,"r_scripts/MD_globular_protein_ligand/ ",part_start,v_list_proteins[i],"/"),ignore.stdout=T,wait = T)
#  system(command = paste0("cp -r ",part_start,"start/toppar/ ",part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/toppar/"),ignore.stdout=T,wait = T)
##  system(command = paste0("cp -r ",part_start,"start/receptor_ligand_MD/resname_mutate.csv ",part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/"),ignore.stdout=T,wait = T)
#  system(command = paste0("cp -r ",part_start,"r_scripts/MD_globular_protein_ligand/ ",part_start,v_list_proteins[i],"/"),ignore.stdout=T,wait = T)
  
}
i<-1
for (i in 1:length(v_list_proteins)) {
  if(file.exists(paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/start/receptor/start.pdb"))){
    part_prepare<-paste0(part_start,",",v_list_proteins[i])
    part_name<-paste0(part_start,v_list_proteins[i],"/MD_globular_protein_ligand/")
    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein_ligand/r_scripts/convert_docking_MD_receptor_ligand.R ",part_prepare),ignore.stdout=T,wait = T)
    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/MD_globular_protein_ligand/r_scripts/prepare_to_stabilisation_MD.R ",part_name),ignore.stdout=T,wait = T)
    #namd
    print(paste0("#run namd2 from "))
    print(paste0(part_name,"namd_script.txt"))
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
