part_analysis <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)
v_rmsds<-seq(from=0,to=100,by=1)
setwd(part_analysis)
part<-paste0(part_analysis,"din/")
setwd(part)

df_all<-read.csv(paste0(part_analysis,"df_all_surf.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand))
if(!dir.exists("df_RMSD_merge_surf")){dir.create("df_RMSD_merge_surf")}
df_all$center<-NULL
df_all<-unique(df_all)
for (i in 1:nrow(df_all)) {
  if(!file.exists(paste0("df_RMSD_merge_surf/",df_all$name[i],".csv"))){
    if(file.exists(paste0("RMSD_merged_surf/",df_all$name[i],".csv"))){
      df_RMSD_all<-read.csv(paste0("RMSD_merged_surf/",df_all$name[i],".csv"),stringsAsFactors = F)
      df_RMSD_all<-df_RMSD_all%>%mutate(RMSD=round(RMSD,digits = 1))
      df_RMSD_all<-df_RMSD_all%>%group_by(RMSD)%>%mutate(number=n())
      df_RMSD_all<-df_RMSD_all%>%select(RMSD,number)
      df_RMSD_all<-unique(df_RMSD_all)
      write.csv(df_RMSD_all,paste0("df_RMSD_merge_surf/",df_all$name[i],".csv"),row.names = F)
    }
  }
}
for (i in 1:nrow(df_all)) {
  if(!file.exists(paste0("df_RMSD_merge_surf/",df_all$name[i],".csv"))){
    df_all$name[i]<-NA
  }
}
df_all<-df_all%>%filter(!is.na(name))
df_RMSD_all<-read.csv(paste0("df_RMSD_merge_surf/",df_all$name[1],".csv"),stringsAsFactors = F)
df_RMSD_all<-df_RMSD_all%>%mutate(name=df_all$name[1])
for (i in 2:nrow(df_all)) {
  df_RMSD_add<-read.csv(paste0("df_RMSD_merge_surf/",df_all$name[i],".csv"),stringsAsFactors = F)
  df_RMSD_add<-df_RMSD_add%>%mutate(name=df_all$name[i])
  df_RMSD_all<-rbind(df_RMSD_all,df_RMSD_add)
}
#df_RMSD_all<-df_RMSD_all%>%group_by(RMSD)%>%mutate(group=sum(number))
#df_RMSD_all<-df_RMSD_all%>%select(RMSD,group)
#df_RMSD_all<-unique(df_RMSD_all)

#df_RMSD_all<-df_RMSD_all%>%filter(RMSD>0)
#seqv<-seq(from=0,to=max(df_RMSD_all$RMSD),by=0.1)



#seqv<-seqv[!seqv%in%df_RMSD_all$RMSD]
#df_RMSD_add<-data.frame(matrix(ncol=2,nrow = length(seqv)))
#colnames(df_RMSD_add)<-colnames(df_RMSD_all)
#df_RMSD_add$RMSD<-seqv
#df_RMSD_add$group<-0
#df_RMSD_all<-rbind(df_RMSD_all,df_RMSD_add)

#df_RMSD_all<-df_RMSD_all%>%mutate(number=RMSD*10)
#df_RMSD_all<-df_RMSD_all%>%mutate(numbera=number%/%1)
#df_RMSD_all<-df_RMSD_all%>%group_by(numbera)%>%mutate(RMSD_new=mean(RMSD))
#df_RMSD_all<-df_RMSD_all%>%group_by(numbera)%>%mutate(group_new=mean(group))
#df_RMSD_all<-ungroup(df_RMSD_all)
#df_RMSD_all<-df_RMSD_all%>%select(RMSD_new,group_new)
#df_RMSD_all<-unique(df_RMSD_all)
p<-ggplot(data=df_RMSD_all)+
#  geom_line(aes(x=RMSD,y=number,colour=name))+
  geom_point(aes(x=RMSD,y=number))+
  theme_bw()+
  facet_grid(name~.)+
  scale_x_continuous(breaks = v_rmsds,labels = v_rmsds)
ggsave(p,filename = paste0("calibration_RMSD_merge_structure_surf.png"), width = 24, height = 15, units = c("cm"), dpi = 200 )
