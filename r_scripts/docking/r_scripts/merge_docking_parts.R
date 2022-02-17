part_name <- commandArgs(trailingOnly=TRUE)
#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-1
group_size<-3
setwd(part_name)
df_all<-read.csv(paste0(part_name,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
group_size<-round(nrow(df_all)/100,digits = 0)
part<-paste0(part_name,"din/")
setwd(part)
if(dir.exists(paste0(part,"fin_merged"))) {system(command = paste0("rm -r ",part,"fin_merged"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"structure_merged"))) {system(command = paste0("rm -r ",part,"structure_merged"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"groups_merged"))) {system(command = paste0("rm -r ",part,"groups_merged"),ignore.stdout=T,wait = T)}
if (!dir.exists("groups_merged")) {dir.create("groups_merged")}
if (!dir.exists("structure_merged")) {dir.create("structure_merged")}
if (!dir.exists("fin_merged")) {dir.create("fin_merged")}

v_protein_name<-unique(df_all$receptor)
v_structure_RMSD<-list.files(paste0("str_fin/"))
df_structure_RMSD<-data.frame(matrix(ncol=4,nrow = length(v_structure_RMSD)))
colnames(df_structure_RMSD)<-c("name","receptor","ligand","RMSD")
df_structure_RMSD$name<-v_structure_RMSD
rm(v_structure_RMSD)
for (j in 1:nrow(df_structure_RMSD)) {
  df_structure_RMSD$receptor[j]<-strsplit(x = df_structure_RMSD$name[j],split = "_")[[1]][1]
  df_structure_RMSD$ligand[j]<-strsplit(x = df_structure_RMSD$name[j],split = "_")[[1]][2]
}
v_ligand<-unique(df_structure_RMSD$ligand)
df_structure_RMSD_TEMP<-df_structure_RMSD%>%filter(ligand==v_ligand[1])
df_structure_RMSD_analysis<-left_join(df_structure_RMSD_TEMP,df_structure_RMSD_TEMP,by=c("receptor","ligand","RMSD"))
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(!is.na(RMSD))
i<-1
j<-1
for (i in 1:length(v_ligand)) {
  df_structure_RMSD_TEMP<-df_structure_RMSD%>%filter(ligand==v_ligand[i])
  df_structure_RMSD_analysis_TEMP<-left_join(df_structure_RMSD_TEMP,df_structure_RMSD_TEMP,by=c("receptor","ligand","RMSD"))
  df_structure_RMSD_analysis_TEMP<-df_structure_RMSD_analysis_TEMP%>%filter(name.x!=name.y)
  
  for (j in 1:nrow(df_structure_RMSD_analysis_TEMP)) {
    pdb_1<-read.pdb(paste0("str_fin/",df_structure_RMSD_analysis_TEMP$name.x[j]))
    pdb_2<-read.pdb(paste0("str_fin/",df_structure_RMSD_analysis_TEMP$name.y[j]))
    
    df_structure_RMSD_analysis_TEMP$RMSD[j]<-rmsd(pdb_1,pdb_2)
  }
  df_structure_RMSD_analysis_TEMP<-df_structure_RMSD_analysis_TEMP%>%filter(RMSD<v_rmsd)
  df_structure_RMSD_analysis<-rbind(df_structure_RMSD_analysis,df_structure_RMSD_analysis_TEMP)
}



df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(RMSD<v_rmsd)
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%group_by(name.x)%>%mutate(number=n())
df_structure_RMSD_analysis<-ungroup(df_structure_RMSD_analysis)       


df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
v_complex<-unique(df_structure_RMSD_analysis$receptor_ligand)
print(v_complex)
i<-1
for (i in 1:length(v_complex)){
  if (!dir.exists(paste0("groups_merged/",v_complex[i]))) {dir.create(paste0("groups_merged/",v_complex[i]))}
  df_structure_RMSD_complex<-df_structure_RMSD_analysis%>%filter(receptor_ligand==v_complex[i])
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
          
          write.csv(df_structure_RMSD_complex_test,paste0("groups_merged/",v_complex[i],"/grop_",j,".csv"),row.names = F) 
        }
      }
      
      df_structure_RMSD_complex$name.x[df_structure_RMSD_complex$name.x%in%c(df_structure_RMSD_complex_test$name.y,df_structure_RMSD_complex_test$name.x)]<-NA
      df_structure_RMSD_complex$name.x[df_structure_RMSD_complex$name.y%in%c(df_structure_RMSD_complex_test$name.y,df_structure_RMSD_complex_test$name.x)]<-NA
      df_structure_RMSD_complex<-df_structure_RMSD_complex%>%filter(!is.na(name.x))
      
      df_RMSD$name.x[df_RMSD$name.x%in%df_structure_RMSD_complex_test$name.y]<-NA
      df_RMSD$name.x[df_RMSD$name.x%in%df_structure_RMSD_complex_test$name.x]<-NA
    }
  }
}

i<-1
v_fin_structure<-list.files("groups_merged/")
for (i in 1:length(v_fin_structure)) {
  v_structure<-list.files(paste0("groups_merged/",v_fin_structure[i]))
  #  if (!dir.exists(paste0("structure_merged/",v_fin_structure[i]))) {dir.create(paste0("structure_merged/",v_fin_structure[i]))}
  df_structure_RMSD_analysis_start<-read.csv(paste0("groups_merged/",v_fin_structure[i],"/",v_structure[1]))
  df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(RMSD<0)
  for (j in 1:length(v_structure)) {
    df_structure_RMSD_analysis<-read.csv(paste0("groups_merged/",v_fin_structure[i],"/",v_structure[j]))
    pdb_1<-read.pdb(paste0("str_fin/",unique(df_structure_RMSD_analysis$name.x)))  
    write.pdb(pdb_1,paste0("structure_merged/",unique(df_structure_RMSD_analysis$name.x)))
    df_structure_RMSD_analysis_start<-rbind(df_structure_RMSD_analysis_start,df_structure_RMSD_analysis)
    
  }
  df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(number>=quantile(df_structure_RMSD_analysis_start$number,probs = 0.25))
  df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%group_by(name.x)%>%mutate(size_of_group=n())
  write.csv(df_structure_RMSD_analysis_start,paste0("fin_merged/",v_fin_structure[i],".csv"),row.names = F)
}

df_structure_RMSD_analysis_start<-read.csv(paste0("fin_merged/",v_fin_structure[1],".csv"),stringsAsFactors = F)
if(length(v_fin_structure)>1){
  for (j in 2:length(v_fin_structure)) {
    df_structure_RMSD_analysis<-read.csv(paste0("fin_merged/",v_fin_structure[j],".csv"),stringsAsFactors = F)
    df_structure_RMSD_analysis_start<-rbind(df_structure_RMSD_analysis_start,df_structure_RMSD_analysis)
  }
}

df_log<-read.csv("log_fin.csv",stringsAsFactors = F)
df_log<-df_log%>%select(models.x,models.y,grop_number,ligand_center,affinity,center,receptor,ligand,new_number)

df_log<-df_log%>%mutate(name=paste0(receptor,"_", ligand,"_",                 center,"_",grop_number,"_",models.x))
df_structure_RMSD<-left_join(df_structure_RMSD_analysis_start,df_log,by=c("name.y" ="name","ligand","receptor"))
df_structure_RMSD<-ungroup(df_structure_RMSD)
df_structure_RMSD<-df_structure_RMSD%>%group_by(name.x)%>%mutate(size_of_group=n())
write.csv(df_structure_RMSD,"df_merge_structure_log.csv",row.names = F)
v_min<-quantile(df_structure_RMSD$affinity,0.025)
v_max<-quantile(df_structure_RMSD$affinity,0.975)
df_structure_RMSD<-df_structure_RMSD%>%group_by(name.x)%>%filter(affinity>v_min)
df_structure_RMSD<-df_structure_RMSD%>%group_by(name.x)%>%filter(affinity<v_max)


a<-seq(from=min(df_structure_RMSD$affinity),to=max(df_structure_RMSD$affinity),by=0.1)
df_structure_RMSD<-df_structure_RMSD%>%mutate(number=as.character(number))
p<-ggplot(data=df_structure_RMSD)+geom_freqpoly(aes(x=affinity,colour=number),binwidth=0.1)+facet_grid(receptor~ligand)+
  scale_x_continuous(breaks=a,labels=a)+theme_bw()+guides(color = "none", size = "none")
ggsave(p,filename = paste0("energy_ligand_receptor_interactions.png"), width = 20, height = 15, units = c("cm"), dpi = 200 )
