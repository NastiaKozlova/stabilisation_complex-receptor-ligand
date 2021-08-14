part_start <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
#part_start<-part_name
part_scriprs<-paste0(part_start,"r_scripts/")
part_start<-paste0(part_start,"docking_first/")
setwd(part_start)

v_ligand<-list.files(paste0("ligand/"))
a<-c()
for (i in 1:length(v_ligand)) {
    b<-strsplit(v_ligand[i],split = ".",fixed = T)[[1]][1]
    a<-c(a,b)
}
v_ligand<-a
  
df_ligand<-data.frame(matrix(ncol=2,nrow=length(v_ligand)))
colnames(df_ligand)<-c("ligand","c")
df_ligand$ligand<-v_ligand
  
df_active_center<-read.csv(paste0("active_center.csv"),stringsAsFactors = F)
v_center<-unique(df_active_center$type)
df_center<-data.frame(matrix(ncol=2,nrow=length(v_center)))
colnames(df_center)<-c("center","c")
df_center$center<-v_center
  
df_ligand_center<-full_join(df_ligand,df_center,by="c")
df_ligand_center$c<-NULL
write.csv(df_ligand_center,paste0("ligand_center.csv"),row.names = F)

#docking
df_ligand_center<-read.csv(paste0("ligand_center.csv"),stringsAsFactors = F)
df_ligand_center<-df_ligand_center%>%mutate(c="C")
v_receptor<-list.files(paste0("receptor_start/"))
a<-c()
for (i in 1:length(v_receptor)) {
  b<-strsplit(v_receptor[i],split = ".",fixed = T)[[1]][1]
  a<-c(a,b)
}
v_receptor<-a
  
df_receptor<-data.frame(matrix(ncol=2,nrow=length(v_receptor)))
colnames(df_receptor)<-c("receptor","c")
df_receptor$receptor<-v_receptor
df_receptor<-df_receptor%>%mutate(c="C")
df_all<-full_join(df_receptor,df_ligand_center,by="c")
df_all$c<-NULL
write.csv(df_all,paste0(part_name,"df_all.csv"),row.names = F)
if (!dir.exists(paste0(part_name,"ligand/"))){dir.create(paste0(part_name,"ligand/"))}
if (!dir.exists(paste0(part_name,"receptor/"))){dir.create(paste0(part_name,"receptor/"))}
i<-1
for (i in 1:nrow(df_receptor)) {
  system(command = paste0(part_start,"programs/MGLTools-1.5.7/bin/pythonsh ",part_start,"programs/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ",
                          part_start,"receptor_start/",df_receptor$receptor[i],".pdb -o ",part_start,"receptor/",df_receptor$receptor[i],".pdbqt ",
                          "-A None"))
}

system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_script.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("chmod +x ",part_start,"script_fin.txt "),ignore.stdout=T,wait = T)
system(command = paste0(part_start,"script_fin.txt"),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"convert_pdbqt_to_pdb.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("chmod +x ",part_start,"convert_pdbqt_to_pdb.py "),ignore.stdout=T,wait = T)
system(command = paste0("python3 ", part_start,"convert_pdbqt_to_pdb.py"),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking/docking_pre_analysis.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking/docking_group_structure.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking/docking_interactions.R ",part_start),ignore.stdout=T,wait = T)
part_start<-paste0(part_start,"")
system(command = paste0("Rscript --vanilla  ",part_start,"docking/analysis.R ",part_start,""),ignore.stdout=T,wait = T)

#setwd(part_name)
#df_topology<-read.csv("din/df_topology.csv",stringsAsFactors = F)
  
#for (i in 1:nrow(df_topology)) {
#  if(!file.exists(paste0("din/groups_fin/",df_topology$name[i],".csv"))){
#    df_topology$name[i]<-NA
#  }
#}
#df_topology<-df_topology%>%filter(!is.na(name))
#df_groups_start<-read.csv(paste0("din/groups_fin/",df_topology$name[1],".csv"),stringsAsFactors = F)
#for (i in 2:nrow(df_topology)) {
#  df_groups_add<-read.csv(paste0("din/groups_fin/",df_topology$name[i],".csv"),stringsAsFactors = F)
#  df_groups_start<-rbind(df_groups_start,df_groups_add)
#}
#write.csv(df_groups_start,"din/df_groups_start.csv",row.names = F)

#df_log<-read.csv(paste0("din/df_log_all.csv"),stringsAsFactors = F)
#df_groups_start<-read.csv("din/df_groups_start.csv",stringsAsFactors = F)
#  colnames(df_groups_start)
#df_groups_start<-df_groups_start%>%select(RMSD, models.y,models.x, number, grop_number, group, ligand_center)
#colnames(df_log)
#df_log<-df_log%>%select(affinity, name_files, name, receptor, ligand, center, files, new_number)
#df_log<-unique(df_log)
#df_log<-df_log%>%mutate(models.y=paste0("frame_",new_number,".pdb"))
#df_logu<-full_join(df_log,df_groups_start,by=c("models.y","name"="ligand_center" ))
#df_logu<-df_logu%>%filter(!is.na(RMSD))
  #  df_logu
#df_logu<-df_logu%>%filter(!is.na(models.y))
#write.csv(df_log,"df_fin_log.csv",row.names = F) 
#df_log<-read.csv("df_fin_log.csv",stringsAsFactors = F)
#df_logu<-df_logu%>%mutate(receptor_fin=NA)
#for (i in 1:nrow(df_logu)) {
#  df_logu$receptor_fin[i]<-strsplit(df_logu$receptor,split = "_")[[1]][1]
#}
#df_logu$number<-NULL
#df_logu<-unique(df_logu)
#p<-ggplot(data=df_logu)+
#  geom_freqpoly(aes(x=affinity))+facet_grid(receptor_fin~ligand)+theme_bw()+guides(colour="none")
#ggsave(p,filename = paste0("din/aminoasids_interated_with_ligands.png"), width = 100, height = 20, units = c("cm"), dpi = 200 ) 
#p<-ggplot(data=df_logu_min)+
#  geom_point(aes(x=number,y=affinity,colour=colour))+facet_grid(receptor~ligand)+theme_bw()+guides(colour="none")
#ggsave(p,filename = paste0("din/aminoasids_interated_with_ligands_low_energy.png"), width = 100, height = 20, units = c("cm"), dpi = 200 ) 

#p<-ggplot(data=df_logu_max)+
#  geom_point(aes(x=number,y=affinity,colour=colour))+facet_grid(receptor~ligand)+theme_bw()+guides(colour="none")
#ggsave(p,filename = paste0("din/aminoasids_interated_with_ligands_higth_energy.png"), width = 100, height = 20, units = c("cm"), dpi = 200 ) 

#df_logu<-df_logu%>%mutate(ligand_center=paste0(ligand,"_",center))
#df_logu<-df_logu%>%mutate(colour="NO")
#df_logu$colour[df_logu$models.y==df_logu$models.x]<-"YES"
#df_logu<-df_logu%>%filter(affinity<0)

#df_logu<-df_logu%>%group_by(name)%>%mutate(min_affinity=min(affinity))
#df_logu<-df_logu%>%group_by(name)%>%mutate(max_affinity=max(affinity))
#df_logu<-df_logu%>%filter(colour=="YES")
#df_logu_min<-df_logu%>%filter(min_affinity==affinity)
#df_logu_max<-df_logu%>%filter(max_affinity==affinity)



#i<-1
#name<-protein[1]
#for (name in protein){
#  part_start_1<-paste0(part, name,"/docking/docking_first/")
#  part_start<-paste0(part,name,"/docking/docking_second/")
#  if(!dir.exists(part_start)){dir.create(part_start)}
#  if(!dir.exists(paste0(part_start,"receptor_start/"))){dir.create(paste0(part_start,"receptor_start/"))}
#  if(!dir.exists(paste0(part_start,"receptor/"))){dir.create(paste0(part_start,"receptor/"))}
  
#  v_start<-list.files(paste0(part_start_1,"din/groups_fin/"))
  #  system(command = paste0("cp -r ",part,name,"/start/ligand_start/ ",part_start,"ligand_start/"),ignore.stdout=T,wait = T)
#  system(command = paste0("cp ",part,name,"/start/active_center.csv ",part_start),ignore.stdout=T,wait = T)
#  df_all<-read.csv(paste0(part_start_1,"df_all.csv"),stringsAsFactors = F)
#  df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
#  for (i in 1:length(v_start)) {
#    df_new_receptor<-read.csv(paste0(part,name,"/docking/docking_first","/din/groups_fin/",v_start[i]),stringsAsFactors = F)
#    #    df_new_receptor<-df_new_receptor%>%filter(number==max(number))
#    df_new_receptor<-df_new_receptor%>%filter(RMSD==0)
#    df_new_receptor<-left_join(df_new_receptor,df_all,by=c("ligand_center"="name"))
#    #    print(paste0(nrow(df_new_receptor)))
#    for (j in 1:nrow(df_new_receptor)) {
#      pdb_lidand<-read.pdb(paste0(part_start_1,"din/pdb_second/",df_new_receptor$ligand_center[j],"/",df_new_receptor$models.x[j]))
#      pdb_receptor<-read.pdb(paste0(part_start_1,"receptor_start/",df_new_receptor$receptor[j],".pdb"))
#      pdb<-cat.pdb(pdb_receptor,pdb_lidand,renumber = T,rechain = T)
#      write.pdb(pdb,paste0(part,name,"/docking/docking_second/receptor_start/",df_new_receptor$ligand_center[j],"_",df_new_receptor$models.x[j]))
#    }
#  }
#  part_start_1<-paste0(part, name,"/docking/docking_first/")
#  part_start<-paste0(part,name,"/docking/docking_second/")
#  df_ligand_center<-read.csv(paste0(part,name,"/start/ligand_center.csv"),stringsAsFactors = F)
#  df_ligand_center<-df_ligand_center%>%mutate(c="C")
#  v_receptor<-list.files(paste0(part_start,"receptor_start/"))
#  a<-c()
#  for (i in 1:length(v_receptor)) {
#    b<-strsplit(v_receptor[i],split = ".",fixed = T)[[1]][1]
#    a<-c(a,b)
#  }
#  v_receptor<-a
#  
#  df_receptor<-data.frame(matrix(ncol=2,nrow=length(v_receptor)))
#  colnames(df_receptor)<-c("receptor","c")
#  df_receptor$receptor<-v_receptor
#  df_receptor<-df_receptor%>%mutate(c="C")
#  df_all<-full_join(df_receptor,df_ligand_center,by="c")
#  df_all$c<-NULL
#  write.csv(df_all,paste0(part_start,"df_all.csv"),row.names = F)
#  if (!dir.exists(paste0(part_start,"ligand/"))){dir.create(paste0(part_start,"ligand/"))}
#  system(command = paste0("cp -r ",part,name,"/start/ligand_start/ ",part_start,"ligand_start/"),ignore.stdout=T,wait = T)
#  system(command = paste0("Rscript --vanilla  ",part,"r_scripts/convert_pdb_to_pdbqt.R ",part_start),ignore.stdout=T,wait = T)
#}

#print("convert pdb to pdbqt receptor-ligand complex in ")

#name<-protein[1]
#for (name in protein){
#  part_start_1<-paste0(part, name,"/docking/docking_first/")
  
#  part_start<-paste0(part,name,"/docking/docking_second/")
#  #  system(command = paste0("Rscript --vanilla  ",part,"r_scripts/docking_script.R ",part_start),ignore.stdout=T,wait = T)
#  #  system(command = paste0("chmod +x ",part_start,"script_fin.txt "),ignore.stdout=T,wait = T)
#  #  print(paste0(part_start,"script_fin.txt"))
#  #  system(command = paste0(part_start,"script_fin.txt"),ignore.stdout=T,wait = T)
#  #  system(command = paste0("Rscript --vanilla  ",part,"r_scripts/convert_pdbqt_to_pdb.R ",part_start),ignore.stdout=T,wait = T)
#  #  system(command = paste0("chmod +x ",part_start,"convert_pdbqt_to_pdb.py "),ignore.stdout=T,wait = T)
#  #  system(command = paste0("python3 ", part_start,"convert_pdbqt_to_pdb.py"),ignore.stdout=T,wait = T)
#  system(command = paste0("Rscript --vanilla  ",part,"r_scripts/docking_pre_analysis.R ",part_start),ignore.stdout=T,wait = T)
#  #to test
#  system(command = paste0("Rscript --vanilla  ",part,"r_scripts/docking_group_structure.R ",part_start),ignore.stdout=T,wait = T)
#  system(command = paste0("Rscript --vanilla  ",part,"r_scripts/docking_interactions.R ",part_start),ignore.stdout=T,wait = T)
#}
#i<-1
#name<-protein[1]
#for (name in protein){
#  #  part_start_1<-paste0(part, name,"/docking/docking_first/")
#  part_start<-paste0(part,name,"/docking/docking_second/")
#  setwd(part_start)
#  df_topology<-read.csv("din/df_topology.csv",stringsAsFactors = F)
  
#  for (i in 1:nrow(df_topology)) {
#    if(!file.exists(paste0("din/groups_fin/",df_topology$name[i],".csv"))){
#      #      print(df_topology$name[i])
#      df_topology$name[i]<-NA
#    }
#  }
#  df_topology<-df_topology%>%filter(!is.na(name))
#  df_groups_start<-read.csv(paste0("din/groups_fin/",df_topology$name[1],".csv"),stringsAsFactors = F)
#  for (i in 2:nrow(df_topology)) {
#    df_groups_add<-read.csv(paste0("din/groups_fin/",df_topology$name[i],".csv"),stringsAsFactors = F)
#    df_groups_start<-rbind(df_groups_start,df_groups_add)
#  }
#  write.csv(df_groups_start,"din/df_groups_start.csv",row.names = F)
#}
#for (name in protein){
#  #  part_start_1<-paste0(part, name,"/docking/docking_first/")
#  part_start<-paste0(part,name,"/docking/docking_second/")
#  setwd(part_start)
#  df_log<-read.csv(paste0(part_start,"din/df_log_all.csv"),stringsAsFactors = F)
#  df_groups_start<-read.csv("din/df_groups_start.csv",stringsAsFactors = F)
#  colnames(df_groups_start)
#  df_groups_start<-df_groups_start%>%select(RMSD, models.y,models.x, number, grop_number, group, ligand_center)
#  colnames(df_log)
#  df_log<-df_log%>%select(affinity, name_files, name, receptor, ligand, center, files, new_number)
#  df_log<-unique(df_log)
#  df_log<-df_log%>%mutate(models.y=paste0("frame_",new_number,".pdb"))
#  df_logu<-full_join(df_log,df_groups_start,by=c("models.y","name"="ligand_center" ))
#  df_logu<-df_logu%>%filter(!is.na(RMSD))
#  #  df_logu
#  df_logu<-df_logu%>%filter(!is.na(models.y))
#  write.csv(df_log,"df_fin_log.csv",row.names = F) 
#  df_logu<-df_logu%>%mutate(ligand_center=paste0(ligand,"_",center))
#  df_logu<-df_logu%>%mutate(colour="NO")
#  df_logu$colour[df_logu$models.y==df_logu$models.x]<-"YES"
#  df_logu<-df_logu%>%filter(affinity<0)
#  
#  df_logu<-df_logu%>%group_by(name)%>%mutate(min_affinity=min(affinity))
#  df_logu<-df_logu%>%group_by(name)%>%mutate(max_affinity=max(affinity))
#  df_logu<-df_logu%>%filter(colour=="YES")
#  df_logu_min<-df_logu%>%filter(min_affinity==affinity)
#  df_logu_max<-df_logu%>%filter(max_affinity==affinity)
#  p<-ggplot(data=df_logu)+
#    geom_point(aes(x=number,y=affinity,colour=colour))+facet_grid(receptor~ligand_center)+theme_bw()+guides(colour="none")
#  ggsave(p,filename = paste0("din/tost.png"), width = 100, height = 40, units = c("cm"), dpi = 200 ) 
#  
  
#  p<-ggplot(data=df_logu_min)+
#    geom_point(aes(x=number,y=affinity,colour=colour))+facet_grid(receptor~ligand_center)+theme_bw()+guides(colour="none")
#  ggsave(p,filename = paste0("din/tost_min.png"), width = 100, height = 40, units = c("cm"), dpi = 200 ) 
  
#  p<-ggplot(data=df_logu_max)+
#    geom_point(aes(x=number,y=affinity,colour=colour))+facet_grid(receptor~ligand_center)+theme_bw()+guides(colour="none")
#  ggsave(p,filename = paste0("din/tost_max.png"), width = 100, height = 40, units = c("cm"), dpi = 200 ) 
#}
#name<-protein[1]
#for (name in protein){
#  #  part_start_1<-paste0(part, name,"/docking/docking_first/")
#  part_start<-paste0(part,name,"/")
#}