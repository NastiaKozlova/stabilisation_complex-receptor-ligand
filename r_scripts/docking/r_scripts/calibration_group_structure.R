part_start <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)
v_rmsds<-seq(from=0,to=100,by=1)
setwd(part_start)
part<-paste0(part_start,"din/")
setwd(part)

df_all<-read.csv(paste0(part_start,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
if(!dir.exists("df_RMSD_all")){dir.create("df_RMSD_all")}
for (i in 1:nrow(df_all)) {
  if(!file.exists(paste0("RMSD_analysis/",df_all$name[i],".csv"))){
    df_all$receptor[i]<-NA
  }
}
df_all<-df_all%>%filter(!is.na(receptor))
for (i in 1:nrow(df_all)) {
  if(!file.exists(paste0("df_RMSD_all/",df_all$name[i],".csv"))){
    df_RMSD_all<-read.csv(paste0("RMSD_analysis/",df_all$name[i],".csv"),stringsAsFactors = F)
    df_RMSD_all<-df_RMSD_all%>%mutate(RMSD=round(RMSD,digits = 1))
    df_RMSD_all<-df_RMSD_all%>%group_by(RMSD)%>%mutate(number=n())
    df_RMSD_all<-df_RMSD_all%>%select(RMSD,number)
    df_RMSD_all<-unique(df_RMSD_all)
    write.csv(df_RMSD_all,paste0("df_RMSD_all/",df_all$name[i],".csv"),row.names = F)
  }
}

df_RMSD_all<-read.csv(paste0("df_RMSD_all/",df_all$name[1],".csv"),stringsAsFactors = F)
for (i in 2:nrow(df_all)) {
  df_RMSD_add<-read.csv(paste0("df_RMSD_all/",df_all$name[i],".csv"),stringsAsFactors = F)
  df_RMSD_all<-rbind(df_RMSD_all,df_RMSD_add)
}
df_RMSD_all<-df_RMSD_all%>%group_by(RMSD)%>%mutate(group=sum(number))
df_RMSD_all<-df_RMSD_all%>%select(RMSD,group)
df_RMSD_all<-unique(df_RMSD_all)

#df_RMSD_all<-df_RMSD_all%>%filter(RMSD>0)
seqv<-seq(from=0,to=max(df_RMSD_all$RMSD),by=0.1)



seqv<-seqv[!seqv%in%df_RMSD_all$RMSD]
df_RMSD_add<-data.frame(matrix(ncol=2,nrow = length(seqv)))
colnames(df_RMSD_add)<-colnames(df_RMSD_all)
df_RMSD_add$RMSD<-seqv
df_RMSD_add$group<-0
df_RMSD_all<-rbind(df_RMSD_all,df_RMSD_add)

df_RMSD_all<-df_RMSD_all%>%mutate(number=RMSD*10)
df_RMSD_all<-df_RMSD_all%>%mutate(numbera=number%/%1)
df_RMSD_all<-df_RMSD_all%>%group_by(numbera)%>%mutate(RMSD_new=mean(RMSD))
df_RMSD_all<-df_RMSD_all%>%group_by(numbera)%>%mutate(group_new=mean(group))
df_RMSD_all<-ungroup(df_RMSD_all)
df_RMSD_all<-df_RMSD_all%>%select(RMSD_new,group_new)
df_RMSD_all<-unique(df_RMSD_all)
p<-ggplot(data=df_RMSD_all)+
  labs(title = "ACHE")+
  geom_line(aes(x=RMSD_new,y=group_new))+
  geom_point(aes(x=RMSD_new,y=group_new))+
  theme_bw()+
  scale_x_continuous(breaks = v_rmsds,labels = v_rmsds)
ggsave(p,filename = paste0("calibration_RMSD_group_structure.png"), width = 24, height = 15, units = c("cm"), dpi = 200 )