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

part_proteins<-paste0(part_start,",",v_list_proteins)
for (i in 1:length(v_list_proteins)) {
  #  part_name<-paste0(part_start,",",v_list_proteins[i])
  part_protein<-part_proteins[i]
  #copying sctipts for docking
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/copy_files_for_docking.R ",part_protein),ignore.stdout=T,wait = T)
  
  #copying sctipts for docking
  #prepare data for docking with and wothout surface counting
  
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/r_scripts/docking_main_surphase_predocking.R ",
                          part_protein),ignore.stdout=T,wait = T)
  #prepare starting files for docking
  part_name<-paste0(part_start,v_list_proteins)[i]
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/r_scripts/prepare_first_docking_main.R ",part_name),ignore.stdout=T,wait = T)
  
}
#}
j<-1
i<-1
part_doking_scriprs<-paste0(part_start,"r_scripts/docking/r_scripts/")
v_search<-c("center","surf")
if(!surphase_conut){
  v_search<-v_search[1]
}
#docking_main_surphase

for (i in 1:length(v_list_proteins)) {
  for (j in 1:length(v_search)) {
    part<-paste0(part_start,v_list_proteins[i],"/docking/docking_first/",v_search[j],"/")
    system(command = paste0("Rscript --vanilla  ",part_doking_scriprs,"docking_script.R ",part),ignore.stdout=T,wait = T)
    
    print("Run")
    print(paste0(part,"script/script_readme.txt"))
  }
}
for (i in 1:length(v_list_proteins)) {
  for (j in 1:length(v_search)) {
    part<-paste0(part_start,v_list_proteins[i],"/docking/docking_first/",v_search[j],"/")
    #    print("Run")
    a<-paste0(part,"script/script_readme.txt")
    if(file.exists(a)){print(a)}
  }
}
i<-1
j<-1
print("Run")
print(paste0(part,"script/",v_receptor[j],"_readme.txt"))
#system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_script.R ",part),ignore.stdout=T,wait = T)
for (i in 1:length(v_list_proteins)) {
  for (j in 1:length(v_search)) {
    
    part_protein<-paste0(part_start,",",v_list_proteins[i],",",v_search[j])
    
    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/r_scripts/docking_main_calculation.R ",part_protein),ignore.stdout=T,wait = T)
    #        system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_script.R ",part),ignore.stdout=T,wait = T)
  }
}
i<-2
for (i in 1:length(v_list_proteins)) {
  for (j in 1:length(v_search)) {
    
    part_protein<-paste0(part_start,",",v_list_proteins[i],",",v_search[j])
    
    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/r_scripts/docking_main_RMSD.R ",part_protein),ignore.stdout=T,wait = T)
    #        system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_script.R ",part),ignore.stdout=T,wait = T)
  }
}
j<-1
i<-1
for (i in 1:length(v_list_proteins)) {
  for (j in 1:length(v_search)) {
    
    part_analysis<-paste0(part_start,v_list_proteins[i],"/docking/docking_first/",v_search[j],"/")
    system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_group_structure.R ",part_analysis,",7"),ignore.stdout=T,wait = T)
    #        system(command = paste0("Rscript --vanilla  ",,"r_scripts/docking/r_scripts/docking_group_structure.R ",
    #                            part_protein,"/docking/","docking_first/",",",1),ignore.stdout=T,wait = T)
  }
}

#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_main_surphase.R ",part_protein),ignore.stdout=T,wait = T)
#}
#}#else{
#  for (i in 1:length(v_list_proteins)) {
#    #prepare castade start_file
#    part_protein<-part_proteins[i]
#    #copying sctipts for docking
#    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_main_center.R ",part_protein),ignore.stdout=T,wait = T)
#    #sctipts for cascade docking
#    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/cascade_docking_center.R ",part_protein),ignore.stdout=T,wait = T)
#  }
#}
