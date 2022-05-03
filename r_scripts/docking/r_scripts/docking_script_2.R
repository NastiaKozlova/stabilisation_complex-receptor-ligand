part_start <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
part_name<-paste0(part_start,"MD_analysis/docking/docking_2/")
setwd(part_name)
if (!dir.exists("tcl")){dir.create("tcl")}
if (!dir.exists("script")){dir.create("script")}
if (!dir.exists("log")){dir.create("log")}
if (!dir.exists("out")){dir.create("out")}
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
df_ligand_center<-read.csv("ligand_center.csv",stringsAsFactors = F)

v_receptor<-unique(df_ligand_center$receptor)
df_active_center_all<-read.csv("active_center.csv",stringsAsFactors = F)
df_active_center<-left_join(df_active_center_all,df_ligand_center,by=c("type"))

j<-1
i<-1
for (j in 1:nrow(df_ligand_center)) {
  pdb<-read.pdb(paste0(part_name,"receptor_start/",df_ligand_center$receptor[j],".pdb"))
  df_pdb<-pdb$atom
  df_pdb<-df_pdb%>%filter(elety=="CA")
  #  start_receptor<-strsplit(v_receptor[j],split = "_")[[1]][1]
  #  df_active_center<-df_active_center_all%>%filter(receptor==start_receptor)
  df_doking<-data.frame(matrix(nrow = 1,ncol = 7))
  colnames(df_doking)<-c('type',"x_len","x_mean","y_len","y_mean","z_len","z_mean")
  df_doking$type<-df_ligand_center$type[j]
  
  for (i in 1:nrow(df_doking)) {
    df_active<-df_active_center%>%filter(type==df_doking$type[i])
    
    df_pdb_a<-df_pdb[df_pdb$resno%in%df_active$resno,]
    df_doking$x_len[i]<-max(df_pdb_a$x)-min(df_pdb_a$x)+10
    df_doking$x_mean[i]<-(max(df_pdb_a$x)+min(df_pdb_a$x))/2
    
    df_doking$y_len[i]<-max(df_pdb_a$y)-min(df_pdb_a$y)+10
    df_doking$y_mean[i]<-(max(df_pdb_a$y)+min(df_pdb_a$y))/2
    
    df_doking$z_len[i]<-max(df_pdb_a$z)-min(df_pdb_a$z)+10
    df_doking$z_mean[i]<-(max(df_pdb_a$z)+min(df_pdb_a$z))/2
  }
  df_doking$x_len[df_doking$x_len>30]<-30
  df_doking$y_len[df_doking$y_len>30]<-30
  df_doking$z_len[df_doking$z_len>30]<-30
  
  #  for (p in 1) {
  for (num in 1:max_num) {
    df_script<-data.frame(matrix(nrow=1,ncol=1))
    df_script[1,1]<-paste0("receptor = receptor/",df_ligand_center$receptor[j],".pdbqt\n",
                           "ligand = ligand/",df_ligand_center$ligand[j],".pdbqt\n\n",
                           "out = out/",df_ligand_center$receptor[j],"_",
                           df_ligand_center$ligand[j],"_",
                           df_doking$type[j],"_",num,".pdbqt\n\n",
                           "center_x = ",df_doking$x_mean[1],
                           "\ncenter_y = ",df_doking$y_mean[1],
                           "\ncenter_z = ",df_doking$z_mean[1],
                           "\n\nsize_x = ", df_doking$x_len[1], 
                           "\nsize_y = ",df_doking$y_len[1], 
                           "\nsize_z = ",df_doking$z_len[1])
    if (!file.exists(paste0("tcl/",df_ligand_center$receptor[j],"_",
                            df_ligand_center$ligand[j],"_",
                            df_doking$type[j],"_",num,"_config.txt"))){
      write.table(df_script,paste0("tcl/",df_ligand_center$receptor[j],"_",
                                   df_ligand_center$ligand[j],"_",
                                   df_ligand_center$type[j],"_",
                                   num,"_config.txt"),row.names = F,col.names = F,sep = "\n",quote = F)
    }
  }
  #}
  df_conf<-data.frame(matrix(nrow = 1,ncol =1))
  a<-paste0(part_start,"programs/autodock_vina_1_1_2_linux_x86/bin/vina --config tcl/",
            df_ligand_center$receptor[j],"_",
            df_ligand_center$ligand[j],"_",
            df_ligand_center$type[j],"_",1:max_num,"_config.txt --log log/",
            df_ligand_center$receptor[j],"_",
            df_ligand_center$ligand[j],"_",
            df_ligand_center$type[j],"_",1:max_num,".log\n")
  a<-paste0(a,collapse = "")
  if (!file.exists(paste0("out/", df_ligand_center$receptor[j],"_",
                          df_ligand_center$ligand[j],"_",
                          df_ligand_center$type[j],"_",
                          max_num,".pdbqt"))){
    df_conf[1,1]<-a
  }
  df_conf_add<-data.frame(matrix(ncol=ncol(df_conf),nrow=1))
  colnames(df_conf_add)<-colnames(df_conf)
  df_conf_add[1,1]<-paste0("cd ",part_name)
  df_conf<-rbind(df_conf_add,df_conf)
  df_conf[is.na(df_conf)]<-""
  write.table(df_conf,paste0("script/",df_ligand_center$receptor[j],"_",
                             df_ligand_center$ligand[j],"_",
                             df_ligand_center$type[j],"_readme.txt"),row.names = F,quote = F,col.names = F,sep = "\n")
  
  system(command = paste0("chmod +x ",part_name,"script/",df_ligand_center$receptor[j],"_",
                          df_ligand_center$ligand[j],"_",df_ligand_center$type[j],"_readme.txt "),ignore.stdout=T,wait = T)
  system(command = paste0(part_name,"script/",df_ligand_center$receptor[j],"_",
                          df_ligand_center$ligand[j],"_",df_ligand_center$type[j],"_readme.txt"),ignore.stdout=T,wait = T)
}
