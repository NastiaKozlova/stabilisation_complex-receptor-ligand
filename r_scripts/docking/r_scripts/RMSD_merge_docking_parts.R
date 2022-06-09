part_analysis <- commandArgs(trailingOnly=TRUE)

#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)

print(Sys.time())

setwd(part_analysis)
df_all<-read.csv(paste0(part_analysis,"df_all.csv"),stringsAsFactors = F)
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
df_all<-df_all%>%filter(!is.na(x))
df_all$name<-NULL

df_active_center<-read.csv("active_center.csv",stringsAsFactors = F)
df_active_center<-left_join(df_all,df_active_center,by=c("center"="type"))

df_analysis<-df_all%>%select(receptor,ligand)
df_analysis<-unique(df_analysis)
df_analysis<-df_analysis%>%mutate(receptor_ligand=paste(receptor,ligand,sep = "_"))
#v_receptor<-unique(df_all$receptor)
j<-1
for (j in 1:nrow(df_analysis)) {
  if(!file.exists(paste0(part_analysis,"din/","RMSD_merged/",df_analysis$receptor_ligand[j],".csv"))){
    df_all_sorted<-df_all%>%filter(receptor==df_analysis$receptor[j])
    df_all_sorted<-df_all_sorted%>%filter(ligand==df_analysis$ligand[j])
    pdb<-read.pdb(paste0(part_analysis,"receptor_start/",df_analysis$receptor[j],".pdb"))
    df_pdb<-pdb$atom
    df_pdb<-df_pdb%>%filter(elety=="CA")
    part<-paste0(part_analysis,"din/")
    setwd(part)
    if (!dir.exists("RMSD_merged")) {dir.create("RMSD_merged")}
    
    df_all_sorted<-df_all_sorted%>%mutate(x_min=NA)
    df_all_sorted<-df_all_sorted%>%mutate(x_max=NA)
    
    df_all_sorted<-df_all_sorted%>%mutate(y_min=NA)
    df_all_sorted<-df_all_sorted%>%mutate(y_max=NA)
    
    df_all_sorted<-df_all_sorted%>%mutate(z_min=NA)
    df_all_sorted<-df_all_sorted%>%mutate(z_max=NA)
    
    for (i in 1:nrow(df_all_sorted)) {
      df_active<-df_active_center%>%filter(center==df_all_sorted$center[i])
      df_pdb_a<-df_pdb[df_pdb$resno%in%df_active$resno,]
      df_all_sorted$x_min[i]<-min(df_pdb_a$x)-5
      df_all_sorted$x_max[i]<-max(df_pdb_a$x)+5
      
      df_all_sorted$y_min[i]<-min(df_pdb_a$y)-5
      df_all_sorted$y_max[i]<-max(df_pdb_a$y)+5
      
      df_all_sorted$z_min[i]<-min(df_pdb_a$z)-5
      df_all_sorted$z_max[i]<-max(df_pdb_a$z)+5
    }
    df_all_sorted_a<-df_all_sorted%>%left_join(df_all_sorted,df_all_sorted,by=c("receptor", "ligand"))
    
    
    df_all_sorted_a<-df_all_sorted_a%>%filter(x.x <=x.y)
    df_all_sorted_a<-df_all_sorted_a%>%filter(y.x <=y.y)
    df_all_sorted_a<-df_all_sorted_a%>%filter(y.x <=z.y)
    
    
    df_all_sorted_a<-df_all_sorted_a%>%filter(x_min.x <=x_min.y)
    df_all_sorted_a<-df_all_sorted_a%>%filter(y_min.x <=y_min.y)
    df_all_sorted_a<-df_all_sorted_a%>%filter(y_min.x <=z_min.y)
    
    df_all_sorted_a<-df_all_sorted_a%>%filter(x_max.x >x_min.y)
    df_all_sorted_a<-df_all_sorted_a%>%filter(y_max.x >y_min.y)
    df_all_sorted_a<-df_all_sorted_a%>%filter(z_max.x >z_min.y)
    df_all_sorted_a<-unique(df_all_sorted_a)
    df_all_sorted_a<-df_all_sorted_a%>%select(receptor,ligand,center.x,center.y)
    
    v_structure<-list.files(paste0(part_analysis,"din/str_fin"))
    df_structure<-data.frame(matrix(nrow=length(v_structure),ncol = 4))
    colnames(df_structure)<-c("name","receptor","ligand","center")
    df_structure$name<-v_structure
    for (i in 1:nrow(df_all_sorted)) {
      df_structure$receptor[grepl(x = df_structure$name,pattern = df_all_sorted$receptor[i])]<-df_all_sorted$receptor[i]
      df_structure$ligand[grepl(x = df_structure$name,pattern = df_all_sorted$ligand[i])]<-df_all_sorted$ligand[i]
      df_structure$center[grepl(x = df_structure$name,pattern = df_all_sorted$center[i])]<-df_all_sorted$center[i]
    }
    
    df_structure<-df_structure%>%filter(!is.na(center))
    #   df_structure<-df_structure
    df_all_sorte<-left_join(df_all_sorted_a,df_structure,by=c("receptor","ligand","center.x"="center"))
    df_all_sorte<-left_join(df_all_sorte,df_structure,by=c("receptor","ligand","center.y"="center"))
    df_all_sorte<-df_all_sorte%>%filter(!is.na(name.x))
    df_all_sorte<-df_all_sorte%>%filter(!is.na(name.y))
    
    df_all_sorte<-df_all_sorte%>%mutate(RMSD=NA)
    df_all_sorte_start<-df_all_sorte%>%filter(!is.na(RMSD))
    print(paste(df_analysis$receptor_ligand[j],Sys.time()))
    v_sorte<-unique(df_all_sorte$center.x)
    for (p in 1:length(v_sorte)) {
      df_all_sorte_add<-df_all_sorte%>%filter(center.x==v_sorte[p])
      df_all_sorte<-df_all_sorte%>%filter(center.x!=v_sorte[p])
      
      for (i in 1:nrow(df_all_sorte_add)) {
        pdb_1<-read.pdb(paste0(part_analysis,"din/str_fin/",df_all_sorte_add$name.x[i]))
        pdb_2<-read.pdb(paste0(part_analysis,"din/str_fin/",df_all_sorte_add$name.y[i]))
        
        df_all_sorte_add$RMSD[i]<-rmsd(pdb_1,pdb_2)
      }
      df_all_sorte_add<-df_all_sorte_add%>%filter(RMSD<50)
      df_all_sorte_start<-rbind(df_all_sorte_start,df_all_sorte_add)
    }
    print(paste(df_analysis$receptor_ligand[j],Sys.time()))
    write.csv(df_all_sorte_start,paste0("RMSD_merged/",df_analysis$receptor_ligand[j],".csv"),row.names=F)
  }
}

print(Sys.time())