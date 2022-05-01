part_name <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-4

setwd(part_name)
setwd("din")

part_name<-paste0(part_start,"MD_analysis/docking/docking_first/")
part<-strsplit(part_name,split = "/")[[1]]
part<-part[1:(length(part)-3)]
part<-paste0(part,collapse = "/")
part<-paste0(part,"/")
df_center<-read.csv(paste0(part,"start/active_center.csv"),stringsAsFactors = F)
df_center<-df_center%>%mutate(center=type)
df_center<-df_center%>%select(center)
df_center<-df_center%>%mutate(c=NA)
df_center<-unique(df_center)
df_all<-read.csv(paste0(part,"start/all_systems.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(receptor=paste0("charmm-gui-",system_name))
df_all<-df_all%>%select(receptor)
df_all<-df_all%>%mutate(c=NA)
df_all<-left_join(df_all,df_center,by="c")
v_ligand<-list.files(paste0(part,"start/ligand_start"))
df_ligand<-data.frame(matrix(nrow = length(v_ligand),ncol=2))
colnames(df_ligand)<-c("ligand","c")
df_ligand$ligand<-v_ligand

for (i in 1:nrow(df_ligand)) {
  df_ligand$ligand[i]<-strsplit(df_ligand$ligand[i],split = ".",fixed = T)[[1]][1]
}
df_all<-left_join(df_all,df_ligand,by="c")
df_all$c<-NULL
df_all<-df_all%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
j<-1
df_all<-df_all%>%mutate(ligand_center=paste0(receptor,"_",ligand,"_",center))

if (dir.exists(paste0("interaction/"))) { system(command = paste0("rm -r ",part_name,"din/interaction/"))}
if (dir.exists(paste0("interaction_TEMP/"))) {system(command = paste0("rm -r ",part_name,"din/interaction_TEMP/"))}
if (dir.exists(paste0("interaction_complex/"))) {system(command = paste0("rm -r ",part_name,"din/interaction_complex/"))}

if (!dir.exists(paste0("interaction/"))) { dir.create(paste0("interaction/"))}
if (!dir.exists("interaction_TEMP/")){dir.create("interaction_TEMP/")}
if (!dir.exists("interaction_complex/")){dir.create("interaction_complex/")}
i<-1
j<-1
p<-1
for (j in 1:nrow(df_all)) {
  if(!file.exists(paste0("groups_fin/",df_all$ligand_center[j],".csv"))){
    df_all$receptor[j]<-NA
  }
}
df_all<-df_all%>%filter(!is.na(receptor))
if(file.exists(paste0("groups_fin/",df_all$ligand_center[j],".csv"))){
  df_model<-read.csv(paste0("groups_fin/",df_all$ligand_center[j],".csv"),stringsAsFactors = F)
  df_model<-df_model%>%group_by(models.x)%>%mutate(size_of_group=n())
}
df_model_all<-read.csv(paste0("groups_fin/",df_all$ligand_center[1],".csv"),stringsAsFactors = F)
df_model_all<-df_model_all%>%mutate(receptor=df_all$receptor[1])
df_model_all<-df_model_all%>%mutate(ligand=df_all$ligand[1])
df_model_all<-df_model_all%>%mutate(center=df_all$center[1])
df_model_all<-df_model_all%>%group_by(models.x)%>%mutate(size_of_group=n())

for (j in 2:nrow(df_all)) {
  
  df_model<-read.csv(paste0("groups_fin/",df_all$ligand_center[j],".csv"),stringsAsFactors = F)
  df_model<-df_model%>%mutate(receptor=df_all$receptor[j])
  df_model<-df_model%>%mutate(ligand=df_all$ligand[j])
  df_model<-df_model%>%mutate(center=df_all$center[j])
  df_model<-df_model%>%group_by(models.x)%>%mutate(size_of_group=n())
  df_model_all<-rbind(df_model_all,df_model)
}
df_all<-df_model_all  
write.csv(df_all,paste0("df_merge_center.csv"),row.names = F)
j<-1
for (j in 1:nrow(df_all)) {
  
  #    df_all<-read.csv(paste0("groups_fin/",df_all$ligand_center[j],".csv"),stringsAsFactors = F)
  #    df_all<-df_all%>%group_by(models.x)%>%mutate(size_of_group=n())
  
  
  a<-read.pdb(paste0(part_name,"receptor_start/",df_all$receptor[j],".pdb"))
  b<-read.pdb(paste0("pdb_second/",df_all$ligand_center[j],"/",df_all$models.y[j]))
  bs<-binding.site(a,b,cutoff=12)
  m<-bs$resnames
  a<-c()
  b<-c()
  y<-1
  for (y in 1:length(m)) {
    p<-strsplit(m[y],split = " ",fixed = T)[[1]][2]
    a<-c(a,p)
    p<-strsplit(m[y],split = " ",fixed = T)[[1]][1]
    b<-c(b,p)
  }
  a<-as.numeric(a)
  df_protein<-data.frame(matrix(ncol=2,nrow=length(a)))
  colnames(df_protein)<-c("resid","resno")
  df_protein$resid<-a
  df_protein$resno<-b
  if (!dir.exists(paste0("interaction/",df_all$ligand_center[j]))) { dir.create(paste0("interaction/",df_all$ligand_center[j]))}
  if (!dir.exists(paste0("interaction/",df_all$ligand_center[j],"/",df_all$size_of_group[j]))) {
    dir.create(paste0("interaction/",df_all$ligand_center[j],"/",df_all$size_of_group[j]))}
  write.csv(df_protein,
            paste0("interaction/",df_all$ligand_center[j],"/",df_all$size_of_group[j],"/",df_all$models.y[j],".csv"),
            row.names = F)
  
}



i<-1
j<-1
df_all<-df_all%>%mutate(receptor_ligand=paste0(receptor,"_",ligand,"_",center))
df_all<-ungroup(df_all)
df_topology<-df_all%>%select(receptor_ligand,receptor,ligand, center,size_of_group)
df_topology<-unique(df_topology)
for (j in 1:nrow(df_topology)) {
  if (!dir.exists(paste0("interaction_TEMP/",df_topology$receptor_ligand[j]))){dir.create(paste0("interaction_TEMP/",df_topology$receptor_ligand[j]))}
  pdb<-read.pdb(paste0(part_name,"receptor_start/",df_topology$receptor[j],".pdb"))
  df_pdb<-pdb$atom
  df_pdb<-df_pdb%>%filter(elety=="CA")
  df_pdb<-df_pdb%>%select(type,resid,resno, x,y,z)
  df_pdb<-df_pdb%>%mutate(number_interactions=0)
  df_pdb<-df_pdb%>%mutate(receptor_ligand=df_topology$receptor_ligand[j])
  df_pdb<-df_pdb%>%mutate(receptor=df_topology$receptor[j])
  df_pdb<-df_pdb%>%mutate(center=df_topology$center[j])
  df_pdb<-df_pdb%>%mutate(ligand=df_topology$ligand[j])
  df_pdb<-df_pdb%>%mutate(size_of_group=df_topology$size_of_group[j])
  v_frame<-list.files(paste0("interaction/",df_topology$receptor_ligand[j],"/",df_topology$size_of_group[j],"/"))
  for (q in 1:length(v_frame)) {
    df_interaction<-read.csv(paste0("interaction/",df_all$receptor_ligand[j],"/",df_all$size_of_group[j],"/",df_all$models.y[j],".csv"),stringsAsFactors = F)
    colnames(df_interaction)<-c(colnames(df_interaction)[2],colnames(df_interaction)[1])
    df_pdb$number_interactions[df_pdb$resno%in%df_interaction$resno]<-df_pdb$number_interactions[df_pdb$resno%in%df_interaction$resno]+1
  }
  df_pdb<-df_pdb%>%mutate(persent_interactions=number_interactions/length(v_frame)*100)
  write.csv(df_pdb,paste0("interaction_TEMP/",df_topology$receptor_ligand[j],"/",df_topology$size_of_group[j],".csv"),row.names = F)
}
i<-1
j<-1

v_system<-list.files(paste0("interaction_TEMP/"))
for (i in 1:length(v_system)) {
  v_frame<-list.files(paste0("interaction_TEMP/",v_system[i]))
  df_pdb<-read.csv(paste0("interaction_TEMP/",v_system[i],"/",v_frame[1]))
  for (q in 2:length(v_frame)) {
    df_pdb_add<-read.csv(paste0("interaction_TEMP/",v_system[i],"/",v_frame[q]))
    df_pdb<-rbind(df_pdb,df_pdb_add)
  }
  df_pdb<-df_pdb%>%mutate(sort=paste(resid,resno,receptor_ligand,size_of_group))
  df_pdb<-df_pdb%>%group_by(sort)%>%mutate(total_interactions=sum(number_interactions))
  df_pdb<-ungroup(df_pdb)
  df_pdb<-df_pdb%>%mutate(total_persent_interactions=total_interactions/size_of_group*100)
  df_pdb<-df_pdb%>%select(resid,resno,
                          receptor_ligand,receptor,ligand,size_of_group,             
                          total_interactions,total_persent_interactions)
  df_pdb<-unique(df_pdb)
  write.csv(df_pdb,paste0("interaction_complex/",v_system[i],".csv"),row.names = F)
}
v_groups<-list.files(paste0("interaction_complex/"))
df_pdb<-read.csv(paste0("interaction_complex/",v_groups[1]),stringsAsFactors =  F)
if(length(v_groups)>1){
  for (i in 2:length(v_groups)) {
    df_pdb_add<-read.csv(paste0("interaction_complex/",v_groups[i]),stringsAsFactors =  F)
    df_pdb<-rbind(df_pdb,df_pdb_add)
  }
}

df_pdb<-df_pdb%>%filter(total_persent_interactions>0)

p<-ggplot(data=df_pdb)+geom_freqpoly(aes(x=total_persent_interactions))+theme_bw()+facet_grid(size_of_group~receptor_ligand)
ggsave(p,filename = paste0(part_name,"interaction_ligand_receptor.png"), width = 12, height = 12, units = c("cm"), dpi = 1000 ) 
