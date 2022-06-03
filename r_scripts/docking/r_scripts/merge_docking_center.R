part_analysis <- commandArgs(trailingOnly=TRUE)
#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)
library(rstatix)
v_rmsd<-1

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

#group_size<-round(nrow(df_all)/100,digits = 0)

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
    v_structure<-c()#unique(df_structure_RMSD_analysis$name.x,df_structure_RMSD_analysis$name.y)
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
    df_RMSD<-df_structure_RMSD_complex%>%select(name.x,number)
    df_RMSD<-unique(df_RMSD)
    df_RMSD<-df_RMSD%>%arrange(desc(number))
    for (j in 1:nrow(df_RMSD)+temp) {
      df_structure_RMSD_complex_test<-df_structure_RMSD_complex_test%>%mutate(grop_number=j)
      df_structure_RMSD_complex_test<-df_structure_RMSD_complex_test%>%filter(grop_number==j)
      write.csv(df_structure_RMSD_complex_test,paste0("groups_merged_center/",df_analysis$receptor_ligand[q],"/grop_",j,".csv"),row.names = F) 
      
    }
  }
  v_structure<-list.files(paste0("groups_merged_center/",df_analysis$receptor_ligand[q]))
  df_structure_RMSD_analysis_start<-read.csv(paste0("groups_merged_center/",df_analysis$receptor_ligand[q],"/",v_structure[1]))
  df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(RMSD<0)
  for (j in 1:length(v_structure)) {
    df_structure_RMSD_analysis<-read.csv(paste0("groups_merged_center/",df_analysis$receptor_ligand[q],"/",v_structure[j]))
    df_structure_RMSD_analysis_start<-rbind(df_structure_RMSD_analysis_start,df_structure_RMSD_analysis)
  }
  
  write.csv(df_structure_RMSD_analysis_start,paste0("fin_merged_center/",df_analysis$receptor_ligand[q],".csv"),row.names = F)
}
df_log<-read.csv("log_fin.csv",stringsAsFactors = F)
df_log<-df_log%>%select(models.x,models.y,grop_number,ligand_center,affinity,center,receptor,ligand,new_number)

df_log<-df_log%>%mutate(name=paste0(receptor,"_", ligand,"_",                 center,"_",grop_number,"_",models.x))
df_log<-df_log%>%mutate(name_add=paste0(receptor,"_", ligand,"_",                 center,"_",grop_number,"_",models.y))
df_log<-df_log%>%group_by(ligand_center)%>%mutate(selection=n())
df_log<-df_log%>%ungroup()
df_log<-df_log%>%select(name,name_add,affinity,selection)
df_structure_RMSD_analysis_start<-read.csv(paste0("fin_merged_center/",df_analysis$receptor_ligand[1],".csv"),stringsAsFactors = F)
#df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%group_by(name.x)%>%mutate(size_of_group=n())
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(is.na(name.x))
j<-1
if(nrow(df_analysis)>0){
  for (j in 1:nrow(df_analysis)) {
    df_structure_RMSD_analysis<-read.csv(paste0("fin_merged_center/",df_analysis$receptor_ligand[j],".csv"),stringsAsFactors = F)
    df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
    df_structure_RMSD_analysis<-rbind(df_structure_RMSD_analysis,df_structure_RMSD_analysis)
    df_structure_RMSD_analysis<-unique(df_structure_RMSD_analysis)
    df_structure_RMSD_analysis_start<-rbind(df_structure_RMSD_analysis_start,df_structure_RMSD_analysis)
  }
}
df_structure_RMSD_analysis_start<-unique(df_structure_RMSD_analysis_start)
df_structure_RMSD_analysis<-left_join(df_structure_RMSD_analysis_start,df_log,by=c("name.y"="name" ))

df_structure_RMSD_analysis<-unique(df_structure_RMSD_analysis)
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%mutate(sort=paste(name.x,center.y))
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%group_by(sort)%>%mutate(center_possibility=n())
df_structure_RMSD_analysis<-ungroup(df_structure_RMSD_analysis)
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%select(name.x,receptor,ligand,center.x,name.y,center.y,       
                                                                          receptor_ligand,name_add,affinity,selection,center_possibility)




df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%group_by(name.x)%>%mutate(group_size=n())
df_structure_RMSD_analysis<-ungroup(df_structure_RMSD_analysis)
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%group_by(name.x)%>%mutate(max_center_possibility=max(center_possibility))
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(max_center_possibility==center_possibility)
df_structure_RMSD_analysis<-ungroup(df_structure_RMSD_analysis)
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%mutate(center=NA)

df_structure_RMSD_analysis<-ungroup(df_structure_RMSD_analysis)

df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%mutate(name=paste(receptor,ligand,center,sep = "_"))
#df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%group_by(name)%>%filter(group_size==max(center_possibility))
df_structure_RMSD_analysis<-ungroup(df_structure_RMSD_analysis)

df_structure<-df_structure_RMSD_analysis%>%filter(center_possibility==max_center_possibility)
df_structure<-df_structure%>%select(name.x,center.y,receptor,ligand,center_possibility,group_size)
df_structure<-unique(df_structure)
for (i in 1:nrow(df_structure)) {
  df_structure_RMSD_analysis$center[df_structure_RMSD_analysis$name.x==df_structure$name.x[i]]<-df_structure$center.y[i]
}


df_structure<-df_structure_RMSD_analysis%>%select(name.x,receptor,ligand,group_size,max_center_possibility,center )
df_structure<-unique(df_structure)
df_structure<-df_structure%>%mutate(name=paste(receptor,ligand,center,sep="_"))
df_structure<-df_structure%>%mutate(sorter=paste(receptor,ligand,center,sep="_"))
df_structure<-ungroup(df_structure)
df_structure<-df_structure%>%group_by(sorter)%>%mutate(max_group_size=max(group_size))
df_structure<-df_structure%>%filter(max_group_size==group_size)
df_structure<-ungroup(df_structure)
length(unique(df_structure$sorter))
for (i in 1:nrow(df_structure)) {
  pdb<-read.pdb(paste0("str_fin/",df_structure$name.x[i]))
  write.pdb(pdb,paste0("structure_merged_center/",df_structure$name[i],"_",df_structure$name.x[i],".pdb"))
}

df_structure_RMSD_analysis<-ungroup(df_structure_RMSD_analysis)

df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%mutate(name=paste(receptor,ligand,center,sep = "_"))

df_structure_RMSD_analysis<-ungroup(df_structure_RMSD_analysis)
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%mutate(group_size=as.character(group_size))

df_structure_RMSD_analysis <- df_structure_RMSD_analysis %>%
  #  filter(center == "cat_trio") %>%
  #  filter(ligand != "PhA") %>%
  mutate(color = ifelse(ligand == "ACh", 1,
                        ifelse(ligand == "choline", 2, 0)))

v_breaks <- c("ACh","ATC","atma","BCh","BTC","BzCh","BzTC", "ontfnac", "PhA2", "Propidium", "TMA", "choline")
v_center<-c("Ацильный карман", "Анионный сайт", "Каталитическая триада", "Оксианионная дыра", "ПАС")
df_center<-df_all%>%select(center)
df_center<-unique(df_center)
df_center<-df_center%>%arrange(center)
df_center<-df_center%>%mutate(center_name=v_center)
df_structure_RMSD_analysis<-left_join(df_structure_RMSD_analysis,df_center,by="center")
v_labels<-c("АХ","АТХ", "АТМА", "БХ", "БТХ", "БзХ", "БзТХ", "о-НФА", "ФА", "Пропидий", "ТМА", "Холин")
#df_structure_RMSD_analysis<-left_join(df_structure_RMSD_analysis,df_center,by="center")
v_test<-seq(from=-10,to=100,by=1)
p<-ggplot(df_structure_RMSD_analysis, aes(x = ligand, y = affinity)) +
  geom_boxplot() +
  labs(y="Аффинность, ККал/моль", x="Лиганд")+
  theme(legend.position = "none") +
  facet_grid(center_name ~ .)+
#  scale_y_continuous(breaks = v_test,labels = v_test) +
#  scale_color_manual(values = c("black", "green", "blue")) +
  theme(legend.position = "none") +
  scale_x_discrete(breaks=v_breaks,labels=v_labels)+
  theme_bw()#+guides(color = "none")
ggsave(p,filename = paste0("energy_ligand_receptor_center.png"), width = 30, height = 20, units = c("cm"), dpi = 200 )
df_structure_RMSD_analysis<-df_structure_RMSD_analysis[df_structure_RMSD_analysis$center%in%c("perefer_anion_cite","cat_trio"),]
p<-ggplot(df_structure_RMSD_analysis, aes(x = ligand, y = affinity)) +
  geom_boxplot() +
  labs(y="Аффинность, ККал/моль", x="Лиганд")+
  theme(legend.position = "none") +
  facet_grid(center_name ~ .)+
  #  scale_y_continuous(breaks = v_test,labels = v_test) +
  #  scale_color_manual(values = c("black", "green", "blue")) +
  theme(legend.position = "none") +
  scale_x_discrete(breaks=v_breaks,labels=v_labels)+
  theme_bw()#+guides(color = "none")
ggsave(p,filename = paste0("energy_ligand_receptor_cat_PAS.png"), width = 20, height = 10, units = c("cm"), dpi = 200 )