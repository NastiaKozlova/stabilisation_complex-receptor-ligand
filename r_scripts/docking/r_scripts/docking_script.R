part <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
setwd(part)
#repetition
max_num<-100

if (!dir.exists("tcl")){dir.create("tcl")}
if (!dir.exists("script")){dir.create("script")}
if (!dir.exists("log")){dir.create("log")}
if (!dir.exists("out")){dir.create("out")}
j<-1
i<-1
df_active_center<-read.csv("active_center.csv",stringsAsFactors = F)
#df_ligand_field<-read.csv("ligand_field.csv",stringsAsFactors = F)
df_all<-read.csv("df_all.csv",stringsAsFactors = F)
j<-1
df_all<-df_all%>%mutate(x_len=NA)%>% mutate(x_mean=NA)%>%
  mutate(y_len=NA)%>%  mutate(y_mean=NA)%>%
  mutate(z_len=NA)%>%  mutate(z_mean=NA)%>%mutate(x_tost=NA)

for (j in 1:nrow(df_all)) {
  pdb<-read.pdb(paste0("receptor_start/",df_all$receptor[j],".pdb"))
  df_pdb<-pdb$atom
  df_pdb<-df_pdb%>%filter(elety=="CA")
  
  df_doking<-df_active_center%>%filter(type==df_all$type[j])
  df_doking<-df_doking
  df_pdb_a<-df_pdb[df_pdb$resno%in%df_doking$resno,]
  
  #df_all$x_tost[j]<-max(df_pdb_a$x)-min(df_pdb_a$x)#+df_all$searching_field[j]
  df_all$x_len[j]<-max(df_pdb_a$x)-min(df_pdb_a$x)+df_all$searching_field[j]
  df_all$x_mean[j]<-(max(df_pdb_a$x)+min(df_pdb_a$x))/2
  
  df_all$y_len[j]<-max(df_pdb_a$y)-min(df_pdb_a$y)+df_all$searching_field[j]
  df_all$y_mean[j]<-(max(df_pdb_a$y)+min(df_pdb_a$y))/2
  
  df_all$z_len[j]<-max(df_pdb_a$z)-min(df_pdb_a$z)+df_all$searching_field[j]
  df_all$z_mean[j]<-(max(df_pdb_a$z)+min(df_pdb_a$z))/2
}
j<-1
for (j in 1:nrow(df_all)) {
  
  for (num in 1:max_num) {
    df_script<-data.frame(matrix(nrow=1,ncol=1))
    df_script[1,1]<-paste0("receptor = receptor/",df_all$receptor[j],".pdbqt\nligand = ligand/",df_all$ligand[j],".pdbqt\n\n",
                           "out = out/",df_all$receptor[j],"_",df_all$ligand[j],"_",df_all$center[j],"_",num,".pdbqt\n\n",
                           "center_x = ",df_all$x_mean[j],"\ncenter_y = ",df_all$y_mean[j],"\ncenter_z = ",df_all$z_mean[j],
                           "\n\nsize_x = ",df_all$x_len[j], "\nsize_y = ",df_all$y_len[j], "\nsize_z = ",df_all$z_len[j])
    write.table(df_script,paste0("tcl/",df_all$receptor[j],"_",df_all$ligand[j],"_",df_all$center[j],"_",num,"_config.txt"),row.names = F,col.names = F,sep = "\n",quote = F)
  }
}

#df_conf<-data.frame(matrix(nrow = length(ligand),ncol =nrow(df_doking)))
df_all<-df_all%>%mutate(script=NA)
for (j in 1:nrow(df_all)) {
  a<-paste0(part,"programs/autodock_vina_1_1_2_linux_x86/bin/vina --config tcl/",
            df_all$receptor[j],"_",df_all$ligand[j],"_",df_all$center[j],"_",1:max_num,"_config.txt --log log/",
            df_all$receptor[j],"_",df_all$ligand[j],"_",df_all$center[j],"_",1:max_num,".log\n")
  a<-paste0(a,collapse = "")
  if ((!file.exists(paste0("out/",df_all$receptor[j],"_",df_all$ligand[j],"_",df_all$center[j],"_",max_num,".pdbqt")))|
      (!file.exists(paste0("log/",df_all$receptor[j],"_",df_all$ligand[j],"_",df_all$center[j],"_",max_num,".log")))){
    df_all$script[j]<-a
  }
}
df_conf<-df_all%>%select(script)
df_conf[is.na(df_conf)]<-""
df_conf_add<-data.frame(matrix(nrow=1,ncol=ncol(df_conf)))
colnames(df_conf_add)<-colnames(df_conf)
df_conf_add[1,1]<-paste0("cd ",part)
df_conf<-rbind(df_conf_add,df_conf)
write.table(df_conf,paste0("script/script_readme.txt"),row.names = F,quote = F,col.names = F,sep = "\n",na="")
system(command = paste0("chmod +x ",part,"script/script_readme.txt"),ignore.stdout=T,wait = T)
#  system(command = paste0(part,"script/",v_receptor[j],"_readme.txt"),ignore.stdout=T,wait = T)
print(paste0(part,"script/script_readme.txt"))

