part_analysis <- commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)
library(rstatix)
setwd(part_analysis)
setwd("din")
df_structure_RMSD<-read.csv("df_merge_structure_log_center.csv",stringsAsFactors =  F)

#df_structure_RMSD<-df_structure_RMSD%>%filter(receptor=="start")
#df_structure_RMSD<-df_structure_RMSD%>%filter(ligand=="ACh")
data_down_results_kw <- df_structure_RMSD %>%
  kruskal_test(affinity ~ name.x)

data_down_results_pwc <- df_structure_RMSD %>%
  tukey_hsd(affinity ~ name.x, p.adjust.method = "BH")
data_down_results_pwc<-ungroup(data_down_results_pwc)
df_structure_RMSD<-df_structure_RMSD%>%mutate(name=paste(receptor,ligand,center,sep = "_"))
df_structure_RMSD<-df_structure_RMSD%>%group_by(name.x)%>%mutate(group_size=n())
df_structure_RMSD<-ungroup(df_structure_RMSD)
df_structure_RMSD<-df_structure_RMSD%>%select(name.x,group_size,receptor,ligand)
df_structure_RMSD<-unique(df_structure_RMSD)
data_down_results_pwc<-left_join(data_down_results_pwc,df_structure_RMSD,by=c("group1"="name.x"))
data_down_results_pwc<-left_join(data_down_results_pwc,df_structure_RMSD,by=c("group2"="name.x"))
data_down_results_pwc<-ungroup(data_down_results_pwc)

data_down_results_pwc_m<-data_down_results_pwc%>%filter(receptor.x==receptor.y)
data_down_results_pwc_m<-data_down_results_pwc_m%>%filter(ligand.x==ligand.y)

data_down_results_pwc_n<-data_down_results_pwc%>%filter(p.adj<0.05)

data_down_results_pwc_pos<-data_down_results_pwc%>%filter(conf.high<0)
data_down_results_pwc_pon<-data_down_results_pwc%>%filter(conf.low>0)
data_down_results_pwc_n<-rbind(data_down_results_pwc_pos,data_down_results_pwc_pon)
