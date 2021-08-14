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
part_name<-paste0(part,"/")
num_model<-3
#part_start<-paste0(part,"/",based_name[name],'MD/stabilisation/')
name<-1
i<-1
for (name in 1:length(based_name)) {
  
  setwd(paste0(part_start,"MD/stabilisation/"))
  #create additional directories
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din'))) {dir.create(paste0(part_start,'MD/stabilisation/din'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/tcl'))) {dir.create(paste0(part_start,'MD/stabilisation/din/tcl'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/pdb_second'))) {dir.create(paste0(part_start,'MD/stabilisation/din/pdb_second'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/pdb_second/',based_name[name]))) {dir.create(paste0(part_start,'MD/stabilisation/din/pdb_second/',based_name[name]))}
  
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/pdb_second_complex'))) {dir.create(paste0(part_start,'MD/stabilisation/din/pdb_second_complex'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/pdb_second_complex/',based_name[name]))) {dir.create(paste0(part_start,'MD/stabilisation/din/pdb_second_complex/',based_name[name]))}
  
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/interaction'))) {dir.create(paste0(part_start,'MD/stabilisation/din/interaction'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/interaction/',based_name[name]))) {dir.create(paste0(part_start,'MD/stabilisation/din/interaction/',based_name[name]))}
  

  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/Energy'))) {dir.create(paste0(part_start,'MD/stabilisation/din/Energy'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/SASA'))) {dir.create(paste0(part_start,'MD/stabilisation/din/SASA'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/RMSD'))) {dir.create(paste0(part_start,'MD/stabilisation/din/RMSD'))}
  if (!dir.exists(paste0(part_start,'MD/stabilisation/din/RMSF'))) {dir.create(paste0(part_start,'MD/stabilisation/din/RMSF'))}
  df_ramachadran<-read.csv(paste0("din/pdb_second/",based_name[name], "_time_Ramachadran.csv"),stringsAsFactors = F)
  df_Energy<-read.table(paste0("din/Energy/interactions_",based_name[name],".txt"), sep="", header=T, na.strings ="", stringsAsFactors= F)
  df_fin<-full_join(df_ramachadran,df_Energy,by=c("number"="Frame"))
  df_fin<-df_fin%>%select(number,ramachadran,Elec,VdW,Nonbond,Total)
  df_fin<-df_fin%>%mutate(RMSD=NA)
  df_full<-full_join(df_fin,df_fin,by="RMSD")
  df_full$time<-NULL
  df_full<-df_full%>%filter(number.x<number.y)
  for (i in 1:nrow(df_full)) {
    pdb_1<-read.pdb(paste0("din/pdb_second_complex/",based_name[name],"/frame_",df_full$number.x[i],".pdb"))
    pdb_2<-read.pdb(paste0("din/pdb_second_complex/",based_name[name],"/frame_",df_full$number.y[i],".pdb"))
    df_full$RMSD[i]<-rmsd(pdb_1,pdb_2,fit = T)
  }
  df_full<-df_full%>%filter(RMSD<5)
  df_full_add<-df_full%>%select(number.y, ramachadran.y, Elec.y, VdW.y, Nonbond.y,  Total.y,  RMSD,   
                                number.x, ramachadran.x, Elec.x, VdW.x, Nonbond.x,  Total.x)
  for (i in 1:nrow(df_fin)) {
    pdb<-read.pdb(paste0("din/pdb_second_complex/",based_name[name],"/frame_",df_fin$number[i],".pdb"))
    a.inds<-atom.select(pdb,"protein")
    b.inds<-atom.select(pdb,"ligand")
    bs<-binding.site(pdb,a.inds = a.inds,b.inds = b.inds)
    m<-bs$resnames
    a<-c()
    b<-c()
    y<-1
    if(length(m)>0){
      for (y in 1:length(m)) {
        p<-strsplit(m[y],split = " ",fixed = T)[[1]][2]
        a<-c(a,p)
        p<-strsplit(m[y],split = " ",fixed = T)[[1]][1]
        b<-c(b,p)
      }
      a<-as.numeric(a)
      df_protein<-data.frame(matrix(ncol=2,nrow=length(a)))
      colnames(df_protein)<-c("resid","resno")
      df_protein$resid<-a
      df_protein$resno<-b
      write.csv(df_protein,paste0("din/interaction/",based_name[name],"/frame_",df_fin$number[i],".csv"),row.names = F)
    }
  }
  v_interactions<-list.files(paste0("din/interaction/",based_name[name],"/"))
  pdb<-read.pdb(paste0("din/pdb_second_complex/",based_name[name],"/frame_",df_fin$number[1],".pdb"))
  df_intecations<-pdb$atom
  df_intecations<-df_intecations%>%mutate(interactions=0)
  df_intecations<-df_intecations%>%filter(elety=="CA")
  i<-1
  for (i in 1:length(v_interactions)) {
    df_intecations_TEMP<-read.csv(paste0("din/interaction/",based_name[name],"/",v_interactions[i]),stringsAsFactors = F)
    df_intecations$interactions[df_intecations$resno%in%df_intecations_TEMP$resid]<-df_intecations$interactions[df_intecations$resno%in%df_intecations_TEMP$resid]+1
  }
  df_intecations<-df_intecations%>%mutate(persent_interactions=interactions/nrow(df_fin)*100)
  write.csv(df_intecations,paste0("din/interaction_",based_name[name],".csv"),row.names = F)
#  v_break<-unique(seq(from=round(min(df_fin$Total),digits = -1),to=round(max(df_fin$Total),digits = -1),by=10))
#  p<-ggplot(data = df_fin)+
#    geom_freqpoly(aes(x=Total),binwidth=1)+theme_bw()+
#    labs(x="Энергия взаимодействия комплекса рецептор-лиганд, ккал/моль",у=" ")+
#    scale_x_continuous(breaks = v_break,labels = v_break)
#  ggsave(p,filename = paste0("energy_interactions_",based_name[name],".png"), width = 20, height = 15, units = c("cm"), dpi = 200 ) 
}

