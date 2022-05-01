part_name <- commandArgs(trailingOnly=TRUE)
#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-3.5

setwd(part_name)
df_all<-read.csv(paste0(part_name,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))

#group_size<-round(nrow(df_all)/100,digits = 0)

part<-paste0(part_name,"din/")
setwd(part)
if(dir.exists(paste0(part,"fin_merged"))) {system(command = paste0("rm -r ",part,"fin_merged"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"structure_merged"))) {system(command = paste0("rm -r ",part,"structure_merged"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"groups_merged"))) {system(command = paste0("rm -r ",part,"groups_merged"),ignore.stdout=T,wait = T)}
if (!dir.exists("RMSD_merged")) {dir.create("RMSD_merged")}
if (!dir.exists("groups_merged")) {dir.create("groups_merged")}
if (!dir.exists("structure_merged")) {dir.create("structure_merged")}
if (!dir.exists("fin_merged")) {dir.create("fin_merged")}

v_protein_name<-unique(df_all$receptor)
v_structure_RMSD<-list.files(paste0("str_fin/"))
df_structure_RMSD<-data.frame(matrix(ncol=4,nrow = length(v_structure_RMSD)))
colnames(df_structure_RMSD)<-c("name","receptor","ligand","RMSD")
df_structure_RMSD$name<-v_structure_RMSD

for (j in 1:nrow(df_structure_RMSD)) {
  df_structure_RMSD$receptor[j]<-strsplit(x = df_structure_RMSD$name[j],split = "_")[[1]][1]
  df_structure_RMSD$ligand[j]<-strsplit(x = df_structure_RMSD$name[j],split = "_")[[1]][2]
}
df_analysis<-df_structure_RMSD%>%select(receptor,ligand)
df_analysis<-unique(df_analysis)
df_analysis<-df_analysis%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
q<-1
for (q in 1:nrow(df_analysis)) {
  df_structure_RMSD_TEMP<-df_structure_RMSD%>%filter(receptor==df_analysis$receptor[q])
  df_structure_RMSD_TEMP<-df_structure_RMSD_TEMP%>%filter(ligand==df_analysis$ligand[q])
  
  group_size<-round(nrow(df_structure_RMSD_TEMP)/100,digits = 0)
  df_structure_RMSD_analysis<-read.csv(paste0("RMSD_merged/",df_analysis$receptor_ligand[q],".csv"),stringsAsFactors = F)
  df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(RMSD<v_rmsd)
  df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%group_by(name.x)%>%mutate(number=n())
  df_structure_RMSD_analysis<-ungroup(df_structure_RMSD_analysis)       
  if (!dir.exists(paste0("groups_merged/",df_analysis$receptor_ligand[q]))) {
    dir.create(paste0("groups_merged/",df_analysis$receptor_ligand[q]))
  }
  
  df_structure_RMSD_complex<-df_structure_RMSD_analysis
  if (nrow(df_structure_RMSD_complex)>0){
    
    df_RMSD<-df_structure_RMSD_complex%>%select(name.x,number)
    df_RMSD<-unique(df_RMSD)
    df_RMSD<-df_RMSD%>%arrange(desc(number))
    df_filter<-df_RMSD
    df_filter<-df_filter%>%mutate(grop_number=1:nrow(df_filter))
    df_structure_RMSD_complex<-df_structure_RMSD_complex
    for (j in 1:nrow(df_RMSD)) {
      if (!is.na(df_RMSD$name.x[j])) {
        df_structure_RMSD_complex_test<-df_structure_RMSD_complex%>%filter(name.x==df_RMSD$name.x[j])
        if (nrow(df_structure_RMSD_complex_test)>group_size) {
          df_structure_RMSD_complex_test<-df_structure_RMSD_complex_test%>%mutate(grop_number=j)
          
          write.csv(df_structure_RMSD_complex_test,paste0("groups_merged/",df_analysis$receptor_ligand[q],"/grop_",j,".csv"),row.names = F) 
        }
      }
      
      df_structure_RMSD_complex$name.x[df_structure_RMSD_complex$name.x%in%c(df_structure_RMSD_complex_test$name.y,df_structure_RMSD_complex_test$name.x)]<-NA
      df_structure_RMSD_complex$name.x[df_structure_RMSD_complex$name.y%in%c(df_structure_RMSD_complex_test$name.y,df_structure_RMSD_complex_test$name.x)]<-NA
      df_structure_RMSD_complex<-df_structure_RMSD_complex%>%filter(!is.na(name.x))
      
      df_RMSD$name.x[df_RMSD$name.x%in%df_structure_RMSD_complex_test$name.y]<-NA
      df_RMSD$name.x[df_RMSD$name.x%in%df_structure_RMSD_complex_test$name.x]<-NA
    }
  }
  v_structure<-list.files(paste0("groups_merged/",df_analysis$receptor_ligand[q]))
  df_structure_RMSD_analysis_start<-read.csv(paste0("groups_merged/",df_analysis$receptor_ligand[q],"/",v_structure[1]))
  df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(RMSD<0)
  for (j in 1:length(v_structure)) {
    df_structure_RMSD_analysis<-read.csv(paste0("groups_merged/",df_analysis$receptor_ligand[q],"/",v_structure[j]))
    df_structure_RMSD_analysis_start<-rbind(df_structure_RMSD_analysis_start,df_structure_RMSD_analysis)
  }
  
  write.csv(df_structure_RMSD_analysis_start,paste0("fin_merged/",df_analysis$receptor_ligand[q],".csv"),row.names = F)
}

df_structure_RMSD_analysis_start<-read.csv(paste0("fin_merged/",df_analysis$receptor_ligand[1],".csv"),stringsAsFactors = F)
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%group_by(name.x)%>%mutate(size_of_group=n())
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(is.na(name.x))
j<-1
if(nrow(df_analysis)>0){
  for (j in 1:nrow(df_analysis)) {
    df_structure_RMSD_analysis<-read.csv(paste0("fin_merged/",df_analysis$receptor_ligand[j],".csv"),stringsAsFactors = F)
    df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
    df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%group_by(name.x)%>%mutate(size_of_group=n())
    df_structure_RMSD_analysis_max<-df_structure_RMSD_analysis%>%group_by(receptor_ligand)%>%filter(size_of_group==max(size_of_group))
    df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(size_of_group>(900/4))
    df_structure_RMSD_analysis<-rbind(df_structure_RMSD_analysis,df_structure_RMSD_analysis_max)
    df_structure_RMSD_analysis<-unique(df_structure_RMSD_analysis)
    df_structure_RMSD_analysis_start<-rbind(df_structure_RMSD_analysis_start,df_structure_RMSD_analysis)
  }
}
df_structure_RMSD_analysis<-df_structure_RMSD_analysis_start%>%select(name.x,receptor,ligand,size_of_group)
df_structure_RMSD_analysis<-unique(df_structure_RMSD_analysis)
for (j in 1:nrow(df_structure_RMSD_analysis)) {
  pdb_1<-read.pdb(paste0("str_fin/",df_structure_RMSD_analysis$name.x[j]))  
  write.pdb(pdb_1,paste0("structure_merged/",df_structure_RMSD_analysis$name.x[j]))
}


df_log<-read.csv("log_fin.csv",stringsAsFactors = F)
df_log<-df_log%>%select(models.x,models.y,grop_number,ligand_center,affinity,center,receptor,ligand,new_number)

df_log<-df_log%>%mutate(name=paste0(receptor,"_", ligand,"_",                 center,"_",grop_number,"_",models.x))
df_structure_RMSD<-left_join(df_structure_RMSD_analysis_start,df_log,by=c("name.y" ="name","ligand","receptor"))
df_structure_RMSD<-ungroup(df_structure_RMSD)
df_structure_RMSD<-df_structure_RMSD%>%group_by(name.x)%>%mutate(size_of_group=n())
write.csv(df_structure_RMSD,"df_merge_structure_log.csv",row.names = F)
a<-seq(from=min(df_structure_RMSD$affinity),to=max(df_structure_RMSD$affinity),by=0.1)
p<-ggplot(data=df_structure_RMSD)+geom_freqpoly(aes(x=affinity,colour=name.x),binwidth=0.1)+facet_grid(receptor~ligand)+
  scale_x_continuous(breaks=a,labels=a)+theme_bw()+guides(color = "none", size = "none")
ggsave(p,filename = paste0("energy_ligand_receptor_interactions.png"), width = 24, height = 15, units = c("cm"), dpi = 200 )
