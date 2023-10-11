part <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
setwd(part)
#repetition
max_num<-100
v_receptor<-list.files("receptor_start")
a<-c()
for (i in 1:length(v_receptor)){
  b<-strsplit(v_receptor[i],split = ".",fixed = T)[[1]][1]
  a<-c(a,b)
}
v_receptor<-a
j<-1
i<-1
df_active_center<-read.csv("active_center.csv",stringsAsFactors = F)

#df_active_center<-df_active_center%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
df_active_center<-df_active_center%>%mutate(x=NA)
df_active_center<-df_active_center%>%mutate(y=NA)
df_active_center<-df_active_center%>%mutate(z=NA)
i<-1
for (i in 1:nrow(df_active_center)) {
  a<-strsplit(df_active_center$type[i],split = "_")[[1]]
  df_active_center$x[i]<-as.numeric(a[3])
  df_active_center$y[i]<-as.numeric(a[5])
  df_active_center$z[i]<-as.numeric(a[7])
}
df_active_center<-df_active_center%>%filter(is.na(x))
df_active_center<-df_active_center%>%select(type,amino,resno)
j<-1
for (j in 1:length(v_receptor)) {
  pdb<-read.pdb(paste0("receptor_start/",v_receptor[j],".pdb"))
  df_pdb<-pdb$atom
  df_pdb<-df_pdb%>%filter(elety=="CA")
  df_doking<-data.frame(matrix(nrow = length(unique(df_active_center$type)),ncol = 7))
  colnames(df_doking)<-c('type',"x_len","x_mean","y_len","y_mean","z_len","z_mean")
  df_doking$type<-unique(df_active_center$type)

  for (i in 1:nrow(df_doking)) {
    df_active<-df_active_center%>%filter(type==df_doking$type[i])
    df_pdb_a<-df_pdb[df_pdb$resno%in%df_active$resno,]
    df_doking$x_len[i]<-max(df_pdb_a$x)-min(df_pdb_a$x)+30
    df_doking$x_mean[i]<-(max(df_pdb_a$x)+min(df_pdb_a$x))/2
  
    df_doking$y_len[i]<-max(df_pdb_a$y)-min(df_pdb_a$y)+30
    df_doking$y_mean[i]<-(max(df_pdb_a$y)+min(df_pdb_a$y))/2
  
    df_doking$z_len[i]<-max(df_pdb_a$z)-min(df_pdb_a$z)+30
    df_doking$z_mean[i]<-(max(df_pdb_a$z)+min(df_pdb_a$z))/2
  }
  a<-list.files("ligand/")
  ligand<-c()
  for (i in 1:length(a)) {
    b<-strsplit(a[i],split = ".",fixed = T)[[1]][1]
    ligand<-c(ligand,b)
  }
  if (!dir.exists("tcl")){dir.create("tcl")}
  if (!dir.exists("script")){dir.create("script")}
  if (!dir.exists("log")){dir.create("log")}
  if (!dir.exists("out")){dir.create("out")}
  for (q in 1:length(ligand)) {
    for (p in 1:nrow(df_doking)) {
      for (num in 1:max_num) {
        df_script<-data.frame(matrix(nrow=1,ncol=1))
        df_script[1,1]<-paste0("receptor = receptor/",v_receptor[j],".pdbqt\nligand = ligand/",ligand[q],".pdbqt\n\nout = out/",v_receptor[j],"_",ligand[q],"_",df_doking$type[p],"_",num,".pdbqt\n\n",
                               "center_x = ",df_doking$x_mean[p],"\ncenter_y = ",df_doking$y_mean[p],"\ncenter_z = ",df_doking$z_mean[p],"\n\nsize_x = ",
                               df_doking$x_len[p], "\nsize_y = ",df_doking$y_len[p], "\nsize_z = ",df_doking$z_len[p])
        if (!file.exists(paste0("tcl/",v_receptor[j],"_",ligand[q],"_",df_doking$type[p],"_",num,"_config.txt"))){
          write.table(df_script,paste0("tcl/",v_receptor[j],"_",ligand[q],"_",df_doking$type[p],"_",num,"_config.txt"),row.names = F,col.names = F,sep = "\n",quote = F)
        }
      }
    }
  }
  df_conf<-data.frame(matrix(nrow = length(ligand),ncol =nrow(df_doking)))
  for (q in 1:length(ligand)) {
    for (p in 1:nrow(df_doking)) {
      a<-paste0(part,"programs/autodock_vina_1_1_2_linux_x86/bin/vina --config tcl/",
                v_receptor[j],"_",ligand[q],"_",df_doking$type[p],"_",1:max_num,"_config.txt --log log/",
                v_receptor[j],"_",ligand[q],"_",df_doking$type[p],"_",1:max_num,".log\n")
      a<-paste0(a,collapse = "")
      if ((!file.exists(paste0("out/", v_receptor[j],"_",ligand[q],"_",df_doking$type[p],"_",max_num,".pdbqt")))|(!file.exists(paste0("log/", v_receptor[j],"_",ligand[q],"_",df_doking$type[p],"_",max_num,".log")))){
        df_conf[q,p]<-a
      }
    }
  }
  df_conf[is.na(df_conf)]<-""
  df_conf_add<-data.frame(matrix(nrow=1,ncol=ncol(df_conf)))
  colnames(df_conf_add)<-colnames(df_conf)
  df_conf_add[1,1]<-paste0("cd ",part)
  df_conf<-rbind(df_conf_add,df_conf)
  write.table(df_conf,paste0("script/",v_receptor[j],"_readme.txt"),row.names = F,quote = F,col.names = F,sep = "\n",na="")
  system(command = paste0("chmod +x ",part,"script/",v_receptor[j],"_readme.txt"),ignore.stdout=T,wait = T)
  system(command = paste0(part,"script/",v_receptor[j],"_readme.txt"),ignore.stdout=T,wait = T)
}
