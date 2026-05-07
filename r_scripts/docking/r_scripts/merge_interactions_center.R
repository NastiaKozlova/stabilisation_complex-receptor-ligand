part_analysis <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-4

setwd(part_analysis)
df_active_center<-read.csv("active_center.csv",stringsAsFactors = F)
df_active_center<-df_active_center%>%mutate(Motif=type)
df_active_center$type<-NULL
setwd("din")


df_all<-read.csv(paste0("log_merge.csv"),stringsAsFactors = F)

df_all<-df_all%>%group_by(ligand,center.x)%>%mutate(min_affinity=min(affinity))
df_all<-df_all%>%filter(min_affinity<0)
p<-ggplot(data=df_all)+
  geom_freqpoly(aes(x=affinity))+
  #  geom_rect(aes(xmin = seq_beg-0.5, xmax = seq_end+0.5, ymin = -Inf, ymax = Inf,fill=topology,alpha=0.1),data=df_topology)+
  facet_grid(ligand~center.x  , scales = "free")+
  #scale_x_continuous(breaks = v_seq,labels = v_seq)+guides(alpha = "none")+
  theme_bw()
p
df_all<-unique(df_all)
df_all<-df_all%>%filter(name.x==name.y)
df_all<-df_all%>%filter(models.x==models.y)
df_all<-unique(df_all)
df_all<-df_all%>%mutate(x_min=NA)
df_all<-df_all%>%mutate(x_max=NA)
df_all<-df_all%>%mutate(y_min=NA)
df_all<-df_all%>%mutate(y_max=NA)
df_all<-df_all%>%mutate(z_min=NA)
df_all<-df_all%>%mutate(z_max=NA)
if (dir.exists(paste0("interaction_merged/"))) {system(command = paste0("rm -r ",part_analysis,"din/interaction_merged/"))}

if (!dir.exists(paste0("interaction_merged/"))) { dir.create(paste0("interaction_merged/"))}
i<-1
j<-1
p<-1

for (j in 1:nrow(df_all)) {
  df_interaction<-read.csv(paste0("interaction/",df_all$ligand_center.y[j],"/",
                                  df_all$new_number[j],".csv"))
  colnames(df_interaction)<-c("resno","resid")
  #df_interaction<-df_interactio
  df_interaction<-df_interaction%>%select(resid,resno)
  df_interaction<-df_interaction%>%mutate(name=df_all$name.x[j])
  df_interaction<-left_join(df_interaction,df_active_center,by=c("resno","resid"="amino"))
  pdb<-read.pdb(paste0("str_fin/",df_all$name.x[j]))
  
  write.pdb(pdb,paste0("str_merged/",df_all$name.x[j]))
  print(paste(unique(df_interaction$name),nrow(df_interaction)))
  write.csv(df_interaction,
            paste0("interaction_merged/",df_all$name.x[j]),row.names = F)
  pdb<-read.pdb(paste0("str_merged/",df_all$name.x[j]))
  df_all$x_min[j]<-min(pdb$atom$x)
  df_all$x_max[j]<-max(pdb$atom$x)
  
  df_all$y_min[j]<-min(pdb$atom$y)
  df_all$y_max[j]<-max(pdb$atom$y)
  
  df_all$z_min[j]<-min(pdb$atom$z)
  df_all$z_max[j]<-max(pdb$atom$z)
}
i<-1

df_interaction<-read.csv(paste0("interaction_merged/",df_all$name.x[i]),stringsAsFactors = F)
if(nrow(df_all)>1){
  for (i in 2:nrow(df_all)) {
    df_interaction_add<-read.csv(paste0("interaction_merged/",df_all$name.x[i]),stringsAsFactors = F)
    df_interaction<-rbind(df_interaction,df_interaction_add)
  }
}


df_all<-df_all%>%select(name.x,center.x,receptor,ligand,x_min,x_max,y_min,y_max,z_min,z_max,affinity)

df_fin<-left_join(df_all,df_interaction,by=c("name.x"="name"))
p<-ggplot(data=df_fin)+
  geom_point(aes(x=center.x,y=Motif))
p
df_HKD<-df_fin[df_fin$Motif%in%c("HKD13","HKD12"),]
df_HKD<-df_HKD%>%filter(ligand=="POPG")
print(unique(df_HKD$name.x))
#df_fin<-left_join(df_fin,df_active_center,by=c("resid"="resno","resno"="amino"),relationship = "many-to-many")
#write.csv(df_pdb,"center_interactions.csv",row.names = F)
#write.csv(df_pdb,"center_interactions.csv",row.names = F)
#p<-ggplot(data=df_fin)+
#  geom_freqpoly(aes(x=affinity))+
  #  geom_rect(aes(xmin = seq_beg-0.5, xmax = seq_end+0.5, ymin = -Inf, ymax = Inf,fill=topology,alpha=0.1),data=df_topology)+
#  facet_grid(ligand~center.x  , scales = "free")+
  #scale_x_continuous(breaks = v_seq,labels = v_seq)+guides(alpha = "none")+
#  theme_bw()
#df_HKD<-df_fin[df_fin$Motif%in%c("HKD13","HKD12"),]

#df_fin_HKD<-df_fin[df_fin$name.x%in%df_HKD$name.x,]
#p<-ggplot(data=df_fin_HKD)+
#  geom_point(aes(x=resid,y=affinity,colour=Motif))+
  #  geom_rect(aes(xmin = seq_beg-0.5, xmax = seq_end+0.5, ymin = -Inf, ymax = Inf,fill=topology,alpha=0.1),data=df_topology)+
#  facet_grid(center.x  ~ligand, scales = "free")+
  #scale_x_continuous(breaks = v_seq,labels = v_seq)+guides(alpha = "none")+
#  theme(axis.text.x = element_blank(), 
#        axis.ticks.x = element_blank()) +
#  theme_bw()
#p
#ggsave(p,   filename = paste0("center_interactions.png"), width = 60, height = 15, units = c("cm"), dpi = 200 ) 
#df_fin_HKD_POPG<-df_fin_HKD%>%filter(ligand=="POPG")
#unique(df_fin_HKD_POPG$name.x)
