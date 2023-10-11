part_analysis <- commandArgs(trailingOnly=TRUE)
name<-strsplit(part_analysis,split = "/",fixed = T)[[1]]
name<-name[length(name)-2]
#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)
library(rstatix)
v_rmsd<-2.5

setwd(part_analysis)
df_all<-read.csv(paste0("df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
df_all<-df_all%>%mutate(x=NA)
df_all<-df_all%>%mutate(y=NA)
df_all<-df_all%>%mutate(z=NA)
i<-1
for (i in 1:nrow(df_all)) {
  a<-strsplit(df_all$center[i],split = "_")[[1]]
  df_all$x[i]<-as.numeric(a[3])
  df_all$y[i]<-as.numeric(a[5])
  df_all$z[i]<-as.numeric(a[7])
}
df_all<-df_all%>%filter(is.na(x))
df_all$name<-NULL

part<-paste0(part_analysis,"din/")
setwd(part)
if(dir.exists(paste0(part,"fin_merged_center"))) {system(command = paste0("rm -r ",part,"fin_merged_center"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"structure_merged_center"))) {system(command = paste0("rm -r ",part,"structure_merged_center"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"groups_merged_center"))) {system(command = paste0("rm -r ",part,"groups_merged_center"),ignore.stdout=T,wait = T)}
if (!dir.exists("RMSD_merged_center")) {dir.create("RMSD_merged_center")}
if (!dir.exists("groups_merged_center")) {dir.create("groups_merged_center")}
if (!dir.exists("structure_merged_center")) {dir.create("structure_merged_center")}
if (!dir.exists("fin_merged_center")) {dir.create("fin_merged_center")}


df_analysis<-df_all%>%select(receptor,ligand)
df_analysis<-unique(df_analysis)
df_analysis<-df_analysis%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
for (q in 1:nrow(df_analysis)) {
  if(!file.exists(paste0("RMSD_merged_center/",df_analysis$receptor_ligand[q],".csv"))){
    df_analysis$receptor[q]<-NA
  }
}
df_analysis<-df_analysis%>%filter(!is.na(receptor))
q<-1
for (q in 1:nrow(df_analysis)) {
  
  df_structure_RMSD_analysis<-read.csv(paste0("RMSD_merged_center/",df_analysis$receptor_ligand[q],".csv"),stringsAsFactors = F)
  df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(RMSD<v_rmsd)
  df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%group_by(name.x)%>%mutate(number=n())
  df_structure_RMSD_analysis<-ungroup(df_structure_RMSD_analysis)       
  if (!dir.exists(paste0("groups_merged_center/",df_analysis$receptor_ligand[q]))) {
    dir.create(paste0("groups_merged_center/",df_analysis$receptor_ligand[q]))
  }
  
  df_structure_RMSD_complex<-df_structure_RMSD_analysis
  if (nrow(df_structure_RMSD_complex)>0){
    
    df_RMSD<-df_structure_RMSD_complex%>%select(name.x,number)
    df_RMSD<-unique(df_RMSD)
    df_RMSD<-df_RMSD%>%arrange(desc(number))
    v_structure<-c()
    df_structure_RMSD_complex<-df_structure_RMSD_complex
    temp<-nrow(df_RMSD)
    for (j in 1:nrow(df_RMSD)) {
      if (!is.na(df_RMSD$name.x[j])) {
        df_structure_RMSD_complex_test<-df_structure_RMSD_complex%>%filter(name.x==df_RMSD$name.x[j])
        if (nrow(df_structure_RMSD_complex_test)>1) {
          df_structure_RMSD_complex_test<-df_structure_RMSD_complex_test%>%mutate(grop_number=j)
          v_structure_TEMP<-unique(c(df_structure_RMSD_complex_test$name.x,df_structure_RMSD_complex_test$name.y))
          v_structure<-c(v_structure,v_structure_TEMP)
          write.csv(df_structure_RMSD_complex_test,paste0("groups_merged_center/",df_analysis$receptor_ligand[q],"/grop_",j,".csv"),row.names = F) 
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
      j<-1
      for (j in (1:nrow(df_RMSD))) {
        df_structure_RMSD_complex_test<-df_structure_RMSD_complex[j,]
        df_structure_RMSD_complex_test<-df_structure_RMSD_complex_test%>%mutate(grop_number=j+temp)
        write.csv(df_structure_RMSD_complex_test,paste0("groups_merged_center/",df_analysis$receptor_ligand[q],"/grop_",j+temp,".csv"),row.names = F) 
        
      }
    }
  }
  v_structure<-list.files(paste0("groups_merged_center/",df_analysis$receptor_ligand[q]))
  df_structure_RMSD_analysis_start<-read.csv(paste0("groups_merged_center/",df_analysis$receptor_ligand[q],"/",v_structure[1]))
  df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(RMSD<0)
  for (j in 1:length(v_structure)) {
    df_structure_RMSD_analysis<-read.csv(paste0("groups_merged_center/",df_analysis$receptor_ligand[q],"/",v_structure[j]))
    df_test<-df_structure_RMSD_analysis%>%filter(is.na(name.x))
    if(nrow(df_test)>0){print(j)}
    df_structure_RMSD_analysis_start<-rbind(df_structure_RMSD_analysis_start,df_structure_RMSD_analysis)
    
  }
  
  write.csv(df_structure_RMSD_analysis_start,paste0("fin_merged_center/",df_analysis$receptor_ligand[q],".csv"),row.names = F)
}
df_structure_RMSD_analysis_start<-read.csv(paste0("fin_merged_center/",df_analysis$receptor_ligand[1],".csv"),stringsAsFactors = F)

df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(is.na(name.x))
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(!is.na(name.x))
j<-1
if(nrow(df_analysis)>0){
  for (j in 1:nrow(df_analysis)) {
    if(file.exists(paste0("fin_merged_center/",df_analysis$receptor_ligand[j],".csv"))){
      df_structure_RMSD_analysis<-read.csv(paste0("fin_merged_center/",df_analysis$receptor_ligand[j],".csv"),stringsAsFactors = F)
      df_structure_RMSD_analysis<-unique(df_structure_RMSD_analysis)
      df_structure_RMSD_analysis_start<-rbind(df_structure_RMSD_analysis_start,df_structure_RMSD_analysis)
    }
  }
}
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%group_by(name.x)%>%mutate(size_of_group=n())
df_structure_RMSD_analysis<-df_structure_RMSD_analysis_start%>%select(name.x,receptor,ligand)

df_log<-read.csv("log_fin_center.csv",stringsAsFactors = F)
df_log<-df_log%>%select(models.x,models.y,grop_number,ligand,affinity,center,receptor,ligand,new_number)

#a<-seq(from=quantile(df_structure_RMSD$affinity,probs = 0.025,na.rm = T),
#       to=quantile(df_structure_RMSD$affinity,na.rm = T,probs = 0.975),by=1)
#a<-round(a,digits = 1)
p<-ggplot(data=df_log)+geom_density(aes(x=affinity,colour=ligand))+facet_grid(center~.)+
  #scale_y_continuous(breaks=a,labels=a)+
  theme_bw()#+guides(color = "none", size = "none")
p
ggsave(p,filename = paste0("energy_ligand_receptor_center.png"), width = 24, height = 15, units = c("cm"), dpi = 200 )


#df_structure_RMSD_analysis<-df_structure_RMSD%>%select(name.x,receptor,ligand,size_of_group)
#df_structure_RMSD_analysis<-unique(df_structure_RMSD_analysis)

df_log<-df_log%>%mutate(name=paste0(receptor,"_", ligand,"_",                 center,"_",grop_number,"_",models.x))
df_structure_RMSD<-left_join(df_structure_RMSD_analysis_start,df_log,by=c("name.x" ="name",
                                                                          "ligand","receptor"),relationship = "many-to-many")
df_structure_RMSD<-ungroup(df_structure_RMSD)
df_structure_RMSD<-df_structure_RMSD%>%group_by(name.x)%>%mutate(size_of_group=n())
write.csv(df_structure_RMSD,"df_merge_structure_log_center.csv",row.names = F)
#a<-seq(from=quantile(df_structure_RMSD$affinity,probs = 0.025,na.rm = T),
#       to=quantile(df_structure_RMSD$affinity,na.rm = T,probs = 0.975),by=1)
#a<-round(a,digits = 1)
#p<-ggplot(data=df_structure_RMSD)+geom_boxplot(aes(y=affinity,x=ligand))+facet_grid(center.x~.)+
#  scale_y_continuous(breaks=a,labels=a)+theme_bw()+guides(color = "none", size = "none")
#ggsave(p,filename = paste0("energy_ligand_receptor_center.png"), width = 24, height = 15, units = c("cm"), dpi = 200 )


#df_structure_RMSD_analysis<-df_structure_RMSD%>%select(name.x,receptor,ligand,size_of_group)
#df_structure_RMSD_analysis<-unique(df_structure_RMSD_analysis)
for (j in 1:nrow(df_structure_RMSD_analysis)) {
  pdb_1<-read.pdb(paste0("str_fin_center/",df_structure_RMSD_analysis$name.x[j]))  
  write.pdb(pdb_1,paste0("structure_merged_center/",df_structure_RMSD_analysis$name.x[j]))
}
