part_start<-"/home/serg/mem/ACHE/"
setwd(part_start)
#install.packages("httr")
#install.packages("cowplot")

#sudo apt-get install libgdal-dev libssl-dev
#sudo apt-get install muscle
#sudo apt-get install openbabel
i<-1
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
#for (i in 1:length(v_list_proteins)) {
#  part_name<-paste0(part_start,v_list_proteins[i],"/docking/")
#  system(command = paste0("Rscript --vanilla  ",part_name,"r_scripts/prepare_first_docking_main.R ",part_name),ignore.stdout=T,wait = T)
#}
i<-1
for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,v_list_proteins[i],"/docking/")
  system(command = paste0("Rscript --vanilla  ",part_name,"r_scripts/first_docking_start_analysis.R ",part_name),ignore.stdout=T,wait = T)
}
for (i in 1:length(v_list_proteins)) {
  part_name<-paste0(part_start,v_list_proteins[i],"/docking/")
  system(command = paste0("Rscript --vanilla  ",part_name,"r_scripts/first_docking_main.R ",part_name),ignore.stdout=T,wait = T)
}
