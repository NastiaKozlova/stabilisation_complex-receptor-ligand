part_name = commandArgs(trailingOnly=TRUE)

# namd run script 
v_namd<-"namd run script"
#quantity of 25 ns MD simulations, change it to Modify length of MD simulation 

num_din<-1
library(dplyr)
library(bio3d)
library(ggplot2)
setwd(part_name)
part<-strsplit(part_name,split = "/",fixed = T)[[1]]
part<-paste0(part[1:(length(part)-1)],collapse = "/")
pdb<-read.pdb(paste0(part,"/docking/docking_first/receptor_start/start.pdb"))
df_pdb<-pdb$atom
df_pdb<-df_pdb%>%filter(elety=="CA")
df_energy<-read.csv("energy_interaction.csv",stringsAsFactors = F)
df_interactions<-read.csv(paste0("prepared_structures/interactions_center_combined.csv"),stringsAsFactors = F)
df_interactions<-df_interactions%>%mutate(name=paste0(pdb_name,".pdb"))
df_combine<-full_join(df_interactions,df_energy,by=c("name"="pdb_structure"))
df_combine<-df_combine%>%mutate(ligand=NA)
for (i in 1:nrow(df_combine)) {
  df_combine$ligand[i]<-strsplit(df_combine$pdb_name[i],split = "_")[[1]][2]
}
df_combine<-left_join(df_combine,df_pdb,by=c("resno","resid"))
df_combine<-df_combine%>%mutate(resid_type="non-polar")
df_combine$resid_type[df_combine$resid%in%c("SER", "THR", "TYR", "ASN", "GLN")]<-"polar"
df_combine$resid_type[df_combine$resid%in%c("LYS", "ARG", "HIS")]<-"positive"
df_combine$resid_type[df_combine$resid%in%c("ASP", "GLU")]<-"negative"
df_combine<-df_combine%>%filter(resid_type!="non-polar")
a<-seq(from=0,to=max(df_combine$resno),by=10)
p_1<-ggplot(data=df_combine)+
  labs(x="Resid number",y="Electrostatic energy, kcal/mol")+
  geom_point(aes(x=resno,y=Elec,colour=resid_type))+
  facet_grid(ligand~.)+
  scale_x_continuous(breaks = a,labels=a)+
  theme(legend.position = "bottom")+
  theme_bw()
#p_1<-ggplot(data=df_arrange)+
#  geom_line(aes(x=mean_amino,y=Elec))+
#  facet_grid(ligand~.)
ggsave(p_1,filename = paste0("energy_interaction.png"), width = 30, height = 20, units = c("cm"), dpi = 200 ) 

df_arrange<-df_combine%>%group_by(pdb_name)%>%mutate(mean_amino=mean(resno))%>%
  mutate(mean_x=mean(x))%>%  mutate(mean_y=mean(y))%>%  mutate(mean_z=mean(z))
  
df_arrange<-ungroup(df_arrange)
df_arrange<-df_arrange%>%select(pdb_name, Elec, VdW, Nonbond, Total,ligand,
                                mean_amino,mean_x,mean_y,mean_z)

df_min<-df_arrange%>%group_by(ligand)%>%mutate(min_Elec=min(Elec))%>%
  mutate(max_Elec=max(Elec))
df_max<-df_min%>%filter(max_Elec==Elec)
df_min<-df_min%>%filter(min_Elec==Elec)

#df_min<-unique(df_min)

p<-ggplot(data=df_arrange)+
  labs(x="Z, A",y="Electrostatic energy, kcal/mol")+
  geom_line(aes(x=mean_z,y=Elec))+
#  geom_smooth(aes(x=mean_z,y=Elec),method="gam")+
  geom_segment(aes(x = mean_z,xend = mean_z,y=-Inf,yend=Inf),colour="red",data = df_min)+
  geom_segment(aes(x = mean_z,xend = mean_z,y=-Inf,yend=Inf),colour="blue",data = df_max)+
#  scale_x_continuous(breaks = seq(from=0,to=max(df_combine$resno,by=10),labels=seq(from=0,to=max(df_combine$resno),by=10))+
  scale_x_continuous(breaks = seq(from=-20,to=20,by=5),labels=seq(from=-20,to=20,by=5))+
  facet_grid(ligand~.)+
  theme_bw()
ggsave(p,filename = paste0("energy_interaction_analysis.png"), width = 30, height = 20, units = c("cm"), dpi = 200 ) 

