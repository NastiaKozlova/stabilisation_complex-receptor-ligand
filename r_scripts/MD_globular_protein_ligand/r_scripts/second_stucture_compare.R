part_start = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
library(ggplot2)
#part_start<-part_analysis
part_name<-part_start
setwd(part_start)
name<-strsplit(part_start,split = "/",fixed = T)[[1]]
based_name<-name[length(name)]
part<-paste0(name[1:(length(name)-1)],collapse = "/")
num_model<-3
#part_start<-paste0(part,"/",based_name[name],'MD/stabilisation/')
name<-1
part_start<-paste0(part,based_name[name],'/MD/stabilisation/')
#find secondaty structure
Extract_Secondary_Structure_From_Pdb<-function(pdb, number_of_pdb){
  m<-dssp(pdb)
  df_second_helix<-data.frame(matrix(nrow = length(m$helix$start),ncol = 3))
  df_second_sheet<-data.frame(matrix(nrow = length(m$sheet$start),ncol = 3))
  if(length(m$sheet$start)>0){
    df_second_sheet$X3<-"sheet"
    for (i in 1:length(m$sheet$start)) {
      df_second_sheet$X1<-m$sheet$start
      df_second_sheet$X2<-m$sheet$end
    }
    df_second_sheet<-df_second_sheet%>%mutate(level_min=number_of_pdb-0.5)
    df_second_sheet<-df_second_sheet%>%mutate(level_max=number_of_pdb+0.5)
  }
  if(length(m$helix$start)>0){
    df_second_helix$X3<-"helix"
    for (i in 1:length(m$helix$start)) {
      df_second_helix$X1<-m$helix$start
      df_second_helix$X2<-m$helix$end
    }
    df_second_helix<-df_second_helix%>%mutate(level_min=number_of_pdb-0.5)
    df_second_helix<-df_second_helix%>%mutate(level_max=number_of_pdb+0.5)
  }
  df_second_all<-rbind(df_second_helix,df_second_sheet)
  colnames(df_second_all)<-c("start","finish","type","level_min","level_max")
  return(df_second_all)
}
j<-1
for (j in 1:length(based_name)) {
  part_start<-paste0(part,"/",based_name[name],'/MD/stabilisation/din/pdb_second')
  setwd(part_start)
  frame_number<-length(list.files(path = paste0(based_name[j])))-1
  pdb<-read.pdb(paste0(based_name[j],"/frame_",0,".pdb"))
  pdb_num<-c((min(pdb$atom$resno)+1):(max(pdb$atom$resno)-1))
  protein.inds <- atom.select(pdb, "protein",resno=pdb_num)
  backpdb <- trim.pdb(pdb, protein.inds)
  df_pdb_all<-Extract_Secondary_Structure_From_Pdb(pdb = backpdb,number_of_pdb = 0)
  for (i in 1:frame_number) {
    pdb_num<-c((min(pdb$atom$resno)+1):(max(pdb$atom$resno)-1))
    protein.inds <- atom.select(pdb, "protein",resno=pdb_num)
    backpdb <- trim.pdb(pdb, protein.inds)
    df_pdb<-Extract_Secondary_Structure_From_Pdb(pdb = backpdb,number_of_pdb = i)
    df_pdb_all<-rbind(df_pdb_all,df_pdb)
    rm(df_pdb)
  }
  write.csv(df_pdb_all,file = paste0("Second_structure_",based_name[j],".csv"),row.names = F)
 #make plot
  p_second<-ggplot(data = df_pdb_all)+
    ggtitle(paste0("Second structure ",based_name[j]))+
    labs(x = "Number of aminoasids", y = "Time (ns)")+
    geom_rect(aes(xmin = start, ymin = level_min, xmax= finish, ymax = level_max, colour = type,fill=type))+
    scale_color_grey()+ scale_fill_grey()+
    theme_bw()+theme(legend.position = "top")
  print(based_name[j])
  ggsave(p_second,filename = paste0("Second_strucure_evolution_",based_name[j],".png"), width = 30, height = 50, units = c("cm"), dpi = 200 ) 
  rm(df_pdb_all)
}
