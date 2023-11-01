part_analysis <- commandArgs(trailingOnly=TRUE)
#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
v_rmsd<-10
v_group_size<-2

setwd(part_analysis)
df_all<-read.csv(paste0(part_analysis,"df_all_surf.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
part<-paste0(part_analysis,"din/")

setwd(part)
if(dir.exists(paste0(part,"fin_merged"))) {system(command = paste0("rm -r ",part,"fin_merged"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"structure_merged_surf"))) {system(command = paste0("rm -r ",part,"structure_merged_surf"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"groups_merged"))) {system(command = paste0("rm -r ",part,"groups_merged"),ignore.stdout=T,wait = T)}
if (!dir.exists("RMSD_merged_surf")) {dir.create("RMSD_merged_surf")}
if (!dir.exists("groups_merged")) {dir.create("groups_merged")}
if (!dir.exists("structure_merged_surf")) {dir.create("structure_merged_surf")}
if (!dir.exists("fin_merged")) {dir.create("fin_merged")}

df_analysis<-df_all%>%select(receptor,ligand)
df_analysis<-unique(df_analysis)
df_analysis<-df_analysis%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
q<-2
for (q in 1:nrow(df_analysis)) {
  if(!file.exists(paste0("RMSD_merged_surf/",df_analysis$receptor_ligand[q],".csv"))){
    df_analysis$receptor_ligand[q]<-NA
  }
}
df_analysis<-df_analysis%>%filter(!is.na(receptor_ligand))
q<-1
for (q in 1:nrow(df_analysis)) {
  
  df_structure_RMSD_analysis<-read.csv(paste0("RMSD_merged_surf/",df_analysis$receptor_ligand[q],".csv"),stringsAsFactors = F)
  df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(RMSD<v_rmsd)
  df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%group_by(name.x)%>%mutate(number=n())
  v_check<-length(unique(c(df_structure_RMSD_analysis$name.x,df_structure_RMSD_analysis$name.y)))
  df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(number>v_check/100)
  df_structure_RMSD_analysis<-ungroup(df_structure_RMSD_analysis)       
  if (!dir.exists(paste0("groups_merged/",df_analysis$receptor_ligand[q]))) {
    dir.create(paste0("groups_merged/",df_analysis$receptor_ligand[q]))
  }
  
  df_structure_RMSD_complex<-df_structure_RMSD_analysis
  if (nrow(df_structure_RMSD_complex)>0){
    
    df_RMSD<-df_structure_RMSD_complex%>%select(name.x,number)
    df_RMSD<-unique(df_RMSD)
    df_RMSD<-df_RMSD%>%arrange(desc(number))
    v_structure<-c()#unique(df_structure_RMSD_analysis$name.x,df_structure_RMSD_analysis$name.y)
    df_structure_RMSD_complex<-df_structure_RMSD_complex
    temp<-nrow(df_RMSD)
    for (j in 1:nrow(df_RMSD)) {
      if (!is.na(df_RMSD$name.x[j])) {
        df_structure_RMSD_complex_test<-df_structure_RMSD_complex%>%filter(name.x==df_RMSD$name.x[j])
        if (nrow(df_structure_RMSD_complex_test)>v_group_size) {
          df_structure_RMSD_complex_test<-df_structure_RMSD_complex_test%>%mutate(grop_number=j)
          v_structure_TEMP<-unique(c(df_structure_RMSD_complex_test$name.x,df_structure_RMSD_complex_test$name.y))
          v_structure<-c(v_structure,v_structure_TEMP)
          write.csv(df_structure_RMSD_complex_test,paste0("groups_merged/",df_analysis$receptor_ligand[q],"/grop_",j,".csv"),row.names = F) 
        }
      }
      
      df_structure_RMSD_complex$name.x[df_structure_RMSD_complex$name.x%in%c(df_structure_RMSD_complex_test$name.y,df_structure_RMSD_complex_test$name.x)]<-NA
      df_structure_RMSD_complex$name.x[df_structure_RMSD_complex$name.y%in%c(df_structure_RMSD_complex_test$name.y,df_structure_RMSD_complex_test$name.x)]<-NA
      df_structure_RMSD_complex<-df_structure_RMSD_complex%>%filter(!is.na(name.x))
      
      df_RMSD$name.x[df_RMSD$name.x%in%df_structure_RMSD_complex_test$name.y]<-NA
      df_RMSD$name.x[df_RMSD$name.x%in%df_structure_RMSD_complex_test$name.x]<-NA
    }
    
    
    df_structure_RMSD_complex<-df_structure_RMSD_analysis
    df_structure_RMSD_complex<-df_structure_RMSD_complex[!df_structure_RMSD_complex$name.x%in%v_structure,]
    df_structure_RMSD_complex<-df_structure_RMSD_complex[!df_structure_RMSD_complex$name.y%in%v_structure,]
    if(nrow(df_structure_RMSD_complex)>0){
      df_RMSD<-df_structure_RMSD_complex%>%select(name.x,number)
      df_RMSD<-unique(df_RMSD)
      df_RMSD<-df_RMSD%>%arrange(desc(number))
      for (j in (1:nrow(df_RMSD))) {
        df_structure_RMSD_complex_test<-df_structure_RMSD_complex[j,]
        df_structure_RMSD_complex_test<-df_structure_RMSD_complex_test%>%mutate(grop_number=j+temp)
        write.csv(df_structure_RMSD_complex_test,paste0("groups_merged/",df_analysis$receptor_ligand[q],"/grop_",j+temp,".csv"),row.names = F) 
      }
    }
  }
  v_structure<-list.files(paste0("groups_merged/",df_analysis$receptor_ligand[q]))
  if(length(v_structure)>0){
    df_structure_RMSD_analysis_start<-read.csv(paste0("groups_merged/",df_analysis$receptor_ligand[q],"/",v_structure[1]))
    df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(RMSD<0)
    for (j in 1:length(v_structure)) {
      df_structure_RMSD_analysis<-read.csv(paste0("groups_merged/",df_analysis$receptor_ligand[q],"/",v_structure[j]))
      df_structure_RMSD_analysis_start<-rbind(df_structure_RMSD_analysis_start,df_structure_RMSD_analysis)
    }
    
    write.csv(df_structure_RMSD_analysis_start,paste0("fin_merged/",df_analysis$receptor_ligand[q],".csv"),row.names = F)
  }
}
if(nrow(df_analysis)>0){
  df_structure_RMSD_analysis_start<-read.csv(paste0("fin_merged/",df_analysis$receptor_ligand[1],".csv"),stringsAsFactors = F)
  df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%group_by(name.x)%>%mutate(size_of_group=n())
  df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
  df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(is.na(name.x))
  j<-1
  if(nrow(df_analysis)>1){
    for (j in 1:nrow(df_analysis)) {
      if(file.exists(paste0("fin_merged/",df_analysis$receptor_ligand[j],".csv"))){
        df_structure_RMSD_analysis<-read.csv(paste0("fin_merged/",df_analysis$receptor_ligand[j],".csv"),stringsAsFactors = F)
        df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
        df_structure_RMSD_analysis<-unique(df_structure_RMSD_analysis)
        df_structure_RMSD_analysis_start<-rbind(df_structure_RMSD_analysis_start,df_structure_RMSD_analysis)
      }
    }
  }
}
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(!is.na(name.x))

df_log<-read.csv("log_fin_surf.csv",stringsAsFactors = F)
df_log<-df_log%>%select(models.x,models.y,grop_number,ligand,affinity,center,receptor,ligand,new_number)

df_log<-df_log%>%mutate(name=paste0(receptor,"_", ligand,"_",                 center,"_",grop_number,"_",models.x))
df_structure_RMSD<-left_join(df_structure_RMSD_analysis_start,df_log,by=c("name.y" ="name","ligand","receptor"))
df_structure_RMSD<-ungroup(df_structure_RMSD)
df_structure_RMSD<-df_structure_RMSD%>%group_by(name.x)%>%mutate(size_of_group=n())
a<-seq(from=0,to=max(df_structure_RMSD$size_of_group),by=20)
p<-ggplot(data = df_structure_RMSD)+
  geom_freqpoly(aes(x=size_of_group),binwidth=1)+
  theme_bw()+
  scale_x_continuous(breaks=a,labels=a)+
  facet_grid(ligand~.)
ggsave(p,filename = paste0("freqpoly_group_size.png"), width = 24, height = 15, units = c("cm"), dpi = 200 )

write.csv(df_structure_RMSD,"df_merge_structure_log_surf.csv",row.names = F)
df_structure_RMSD<-read.csv("df_merge_structure_log_surf.csv",stringsAsFactors = F)

#df_structure_RMSD<-df_structure_RMSD%>%filter(size_of_group>5)
a<-seq(from=round(min(df_structure_RMSD$affinity,na.rm = T),digits = 0),to=round(max(df_structure_RMSD$affinity,na.rm = T),digits = 0),by=0.5)
#df_structure_RMSD<-df_structure_RMSD[df_structure_RMSD$ligand%in%c("DMPE","DPPE", "DYPE","PYPE"),]
df_structure_RMSD_analysis<-ungroup(df_structure_RMSD)
df_log<-df_log%>%group_by(ligand)%>%mutate(max_affinity=quantile(affinity,probs = 0.975))%>%
  mutate(min_affinity=quantile(affinity,probs = 0.025))
df_log<-df_log%>%filter(affinity>min_affinity)%>%filter(affinity<max_affinity)

p<-ggplot(data= df_log)+
  geom_density(aes(x=affinity,colour=ligand))+
  theme_bw()
ggsave(p,   filename = paste0("energy_ligand_receptor_interactions_surf.png"), width = 32, height = 18, units = c("cm"), dpi = 200 ) 
#energy_ligand_receptor_interactions
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(name.x==name.y)
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%select(name.x,receptor,ligand,size_of_group)

df_structure_RMSD_analysis<-unique(df_structure_RMSD_analysis)
for (j in 1:nrow(df_structure_RMSD_analysis)) {
  pdb_1<-read.pdb(paste0("str_fin/",df_structure_RMSD_analysis$name.x[j]))  
  write.pdb(pdb_1,paste0("structure_merged_surf/",df_structure_RMSD_analysis$name.x[j]))
}
