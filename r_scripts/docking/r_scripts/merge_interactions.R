part_name <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-4

setwd(part_name)
setwd("din")


df_all<-read.csv(paste0("df_merge_structure_log.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(receptor_ligand=paste0(receptor,"_",ligand,"_",center))

#if (dir.exists(paste0("interaction/"))) { system(command = paste0("rm -r ",part_name,"din/interaction/"))}
#if (dir.exists(paste0("interaction_TEMP/"))) {system(command = paste0("rm -r ",part_name,"din/interaction_TEMP/"))}
if (dir.exists(paste0("interaction_serf/"))) {system(command = paste0("rm -r ",part_name,"din/interaction_serf/"))}

if (!dir.exists(paste0("interaction_serf/"))) { dir.create(paste0("interaction_serf/"))}
i<-1
j<-3
p<-1
v_structure<-unique(df_all$name.x)
for (j in 1:length(v_structure)) {
  df_complex<-df_all%>%filter(name.x==v_structure[j])
  pdb<-read.pdb(paste0(part_name,"receptor_start/",df_all$receptor[j],".pdb"))
  
  df_pdb<-pdb$atom
  df_pdb<-df_pdb%>%filter(elety=="CA")
  df_pdb<-df_pdb%>%mutate(number_interactions=0)
  df_pdb<-df_pdb%>%mutate(tested_structure=0)
  df_pdb<-df_pdb%>%mutate(total_structure=nrow(df_complex))
  test<-nrow(df_pdb)
  for (p in 1:nrow(df_complex)) {
    if(file.exists(paste0("interaction/",df_complex$receptor_ligand[p],"/",df_complex$new_number[p],".csv"))){
      df_protein<-read.csv(paste0("interaction/",df_complex$receptor_ligand[p],"/",df_complex$new_number[p],".csv"),
                           stringsAsFactors = F) 
      df_pdb$number_interactions[df_pdb$resno%in%df_protein$resid]<-df_pdb$number_interactions[df_pdb$resno%in%df_protein$resid]+1
      df_pdb$tested_structure<-df_pdb$tested_structure+1
      
    }
  }
  df_pdb<-df_pdb%>%filter(tested_structure==total_structure)
  if(nrow(df_pdb)==test){
    df_pdb<-df_pdb%>%select(resno,resid,x,y,z,number_interactions,tested_structure,total_structure)
    df_pdb<-df_pdb%>%mutate(persent_interactions=number_interactions/total_structure*100)
    write.csv(df_pdb,
              paste0("interaction_serf/",v_structure[j],".csv"),row.names = F)
  }
}
df_all<-df_all%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
df_all<-df_all%>%select(name.x,receptor,ligand,receptor_ligand)
df_all<-unique(df_all)
df_all<-df_all%>%group_by(receptor_ligand)%>%mutate(number=1:n())
df_all<-df_all%>%mutate(number=as.character(number))
df_all<-ungroup(df_all)
df_pdb<-read.csv(paste0("interaction_serf/",df_all$name.x[1],".csv"),stringsAsFactors = F)
df_pdb<-df_pdb%>%mutate(name.x=df_all$name.x[1])
df_pdb<-df_pdb%>%filter(persent_interactions==100)
for (j in 2:nrow(df_all)) {
  df_pdb_add<-read.csv(paste0("interaction_serf/",df_all$name.x[j],".csv"),stringsAsFactors = F)
  df_pdb_add<-df_pdb_add%>%mutate(name.x=df_all$name.x[j])
  df_pdb<-rbind(df_pdb,df_pdb_add)
  df_pdb<-df_pdb%>%filter(persent_interactions==100)
}
#df_pdb<-df_pdb%>%filter(persent_interactions==100)
df_pdb<-left_join(df_pdb,df_all,by="name.x")
#df_pdb<-df_pdb[df_pdb$ligand%in%c("Na","HPO4"),]
#df_pdb<-df_pdb%>%filter(receptor=="charmm-gui-1717818438")
part_start<-strsplit(part_name,split = "/",fixed = T)[[1]]
part_start<-paste0(part_start[1:(length(part_start)-3)],collapse = "/")
part_start<-paste0(part_start,"/")
#df_topology<-read.csv(paste0(part_start,"start/df_topology.csv"),stringsAsFactors = F)
#v_seq<-seq(from=0,to =max(df_topology$seq_end),by=10)
v_seq<-seq(from=0,to =10000,by=10)
p<-ggplot()+
#  geom_rect(aes(xmin = seq_beg-0.5, xmax = seq_end+0.5, ymin = -Inf, ymax = Inf,fill=topology,alpha=0.1),data=df_topology)+
  geom_point(aes(x=resno,y=number),data=df_pdb)+
  
  theme_bw()+facet_grid(receptor_ligand~., scales = "free")+
  scale_x_continuous(breaks = v_seq,labels = v_seq)+
  guides(alpha = "none")
ggsave(p,   filename = paste0("serf_interactions.png"), width = 60, height = 30, units = c("cm"), dpi = 200 ) 