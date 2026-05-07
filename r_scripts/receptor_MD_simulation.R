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