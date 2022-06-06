part_analysis <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)
v_rmsds<-seq(from=0,to=100,by=1)
setwd(part_analysis)
part<-paste0(part_analysis,"din/")
setwd(part)

df_all<-read.csv(paste0(part_analysis,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
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
#df_all$name<-NULL
df_all<-df_all%>%select(name,receptor,ligand,center)

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
    df_RMSD_all<-df_RMSD_all%>%mutate(RMSD=round(RMSD,digits = 2))
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

#df_RMSD_all<-ungroup(df_RMSD_all)
#df_RMSD_all<-df_RMSD_all%>%select(RMSD_new,group_new)
#df_RMSD_all<-unique(df_RMSD_all)
p<-ggplot(data=df_RMSD_all)+
#  labs(title = "ACHE")+
  geom_line(aes(x=RMSD,y=group))+
  geom_point(aes(x=RMSD,y=group))+
  theme_bw()+
  scale_x_continuous(breaks = v_rmsds,labels = v_rmsds)
ggsave(p,filename = paste0("calibration_RMSD_group_structure_center.png"), width = 24, height = 15, units = c("cm"), dpi = 200 )
