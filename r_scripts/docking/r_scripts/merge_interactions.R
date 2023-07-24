part_name <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-4

setwd(part_name)
setwd("din")

part_temp<-strsplit(part_name,split = "/",fixed = T)[[1]]
a<-paste0(part_temp[1:(length(part_temp)-3)],collapse = "/")
b<-paste0("/start/sequence/",part_temp[length(part_temp)-2],".fasta",collapse = "")
seq_name<-paste0(a,b)
seq<-read.fasta(seq_name)

v_seq<-as.vector(seq$ali)
df_seq<-data.frame(matrix(nrow = length(v_seq),ncol=2))
colnames(df_seq)<-c("resno","resid")
df_seq$resno<-1:nrow(df_seq)
df_seq$resid<-v_seq
df_all<-read.csv(paste0("df_merge_structure_log.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(receptor_ligand=paste0(receptor,"_",ligand,"_",center))


if (dir.exists(paste0("interaction_surf/"))) {system(command = paste0("rm -r ",part_name,"din/interaction_surf/"))}

if (!dir.exists(paste0("interaction_surf/"))) { dir.create(paste0("interaction_surf/"))}
i<-1
j<-3
p<-1
v_structure<-unique(df_all$name.x)
v_structure_test<-c()
for (j in 1:length(v_structure)) {
   if(file.exists(paste0("structure_merged/",v_structure[j]))){
      v_structure_test<-c(v_structure_test,v_structure[j])
   }
}
df_all<-df_all[df_all$name.x%in%v_structure_test,]
v_structure<-v_structure_test
for (j in 1:length(v_structure)) {
   df_complex<-df_all%>%filter(name.x==v_structure[j])
   pdb<-read.pdb(paste0(part_name,"receptor_start/",df_all$receptor[j],".pdb"))
   
   df_pdb<-pdb$atom
   df_pdb<-df_pdb%>%filter(elety=="CA")
   df_pdb<-df_pdb%>%mutate(number_interactions=0)
   df_pdb<-df_pdb%>%mutate(tested_structure=0)
   df_pdb<-df_pdb%>%mutate(total_structure=nrow(df_complex))
   test<-nrow(df_pdb)
   for (p in 1:nrow(df_complex)) {
      
      if(file.exists(paste0("interaction/",df_complex$receptor_ligand[p],"/",df_complex$new_number[p],".csv"))){
         df_protein<-read.csv(paste0("interaction/",df_complex$receptor_ligand[p],"/",df_complex$new_number[p],".csv"),
                              stringsAsFactors = F) 
         df_pdb$number_interactions[df_pdb$resno%in%df_protein$resid]<-df_pdb$number_interactions[df_pdb$resno%in%df_protein$resid]+1
         df_pdb$tested_structure<-df_pdb$tested_structure+1
         
      }
   }
   df_pdb<-df_pdb%>%filter(tested_structure==total_structure)
   if(nrow(df_pdb)==test){
      df_pdb<-df_pdb%>%select(resno,resid,x,y,z,number_interactions,tested_structure,total_structure)
      df_pdb<-df_pdb%>%mutate(persent_interactions=number_interactions/total_structure*100)
      
      write.csv(df_pdb,
                paste0("interaction_surf/",v_structure[j],".csv"),row.names = F)
   }
}

#df_all_test<-df_all
#df_all<-df_all_test
#df_all<-df_all[df_all$name.x%in%v_structure_test,]
df_all<-df_all%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
df_all<-df_all%>%select(name.x,receptor,ligand,receptor_ligand)
df_all<-unique(df_all)
df_all<-df_all%>%group_by(receptor_ligand)%>%mutate(number=1:n())
df_all<-df_all%>%mutate(number=as.character(number))
df_all<-ungroup(df_all)
df_pdb<-read.csv(paste0("interaction_surf/",df_all$name.x[1],".csv"),stringsAsFactors = F)
df_pdb<-df_pdb%>%mutate(name.x=df_all$name.x[1])
df_pdb<-df_pdb%>%filter(persent_interactions==100)
for (j in 2:nrow(df_all)) {
   df_pdb_add<-read.csv(paste0("interaction_surf/",df_all$name.x[j],".csv"),stringsAsFactors = F)
   df_pdb_add<-df_pdb_add%>%mutate(name.x=df_all$name.x[j])
   df_pdb<-rbind(df_pdb,df_pdb_add)
   df_pdb<-df_pdb%>%filter(persent_interactions==100)
}

df_pdb<-left_join(df_pdb,df_all,by="name.x")

#part_start<-strsplit(part_name,split = "/",fixed = T)[[1]]
#part_start<-paste0(part_start[1:(length(part_start)-3)],collapse = "/")
#part_start<-paste0(part_start,"/")
#df_topology<-read.csv(paste0(part_start,"start/df_topology.csv"),stringsAsFactors = F)
v_seq<-seq(from=0,to =max(df_pdb$resno),by=30)
#v_seq<-seq(from=0,to =10000,by=10)
#"MET" "SER" "LEU" "TRP" "LYS" "ILE" "GLY" "VAL" "ALA" "PHE" "THR" "HIS" "ASP" "ARG" "PRO" "TYR" "GLU" "GLN" "ASN" "CYS"
#df_pdb<-df_pdb[df_pdb$ligand%in%c("DMPE", "DPPG", "DYPE", "DYPG",  "PYPE"),]
#df_pdb<-df_pdb[df_pdb$ligand%in%c("DMPE", "DPPG", "DYPE", "DYPG"),]
df_pdb<-df_pdb%>%group_by(name.x)%>%mutate(total_structure=median(resno))
df_structure<-df_pdb%>%select(name.x, total_structure, receptor,ligand,receptor_ligand)
df_structure<-unique(df_structure)
df_structure<-df_structure%>%group_by(ligand,receptor)%>%arrange(total_structure)%>%mutate(number=1:n())
df_pdb$number<-NULL
df_pdb<-left_join(df_pdb,df_structure)
#df_pdb<-df_pdb%>%group_by(name.x)%>%mutate(total_structure=median(resno))
df_pdb<-df_pdb%>%mutate(number=as.numeric(number))
#write.csv(df_pdb,"full_aminoacids_interactions.csv",row.names = F)
#ggplot(df_pdb)+
#   geom_point(aes(x =total_structure , y =number))
df_pdb$resid[df_pdb$resid%in%c("LYS", "ARG", "HIS")]<-"POSITIVE"
df_pdb$resid[df_pdb$resid%in%c("ASP", "GLU")]<-"NEGATIVE"
df_pdb$resid[df_pdb$resid%in%c("MET","SER", "LEU","ILE", "GLY", "VAL", "ALA", "THR", "PRO","GLN", "ASN", "CYS")]<-"NEUTRAL"
df_pdb$resid[df_pdb$resid%in%c("TRP", "PHE", "TYR")]<-"AROMATIC"
write.csv(df_pdb,"full_aminoacids_interactions.csv",row.names = F)
p<-ggplot(df_pdb)+
   geom_rect(aes(xmin = resno-0.5, xmax = resno+0.5, ymin =number-0.5 , ymax = number+0.5,fill=resid,colour=resid))+
   #geom_freqpoly(aes(x=resno),binwidth=1,data=df_pdb)+
   #geom_point(aes(x =resno , y =number))+
   
   theme_bw()+facet_grid(ligand~receptor, scales = "free")+
   scale_x_continuous(breaks = v_seq,labels = v_seq)+
   guides(alpha = "none")
ggsave(p,   filename = paste0("surf_interactions.png"), width = 10*length(unique(df_pdb$ligand)/12*9, height = 10*length(unique(df_pdb$ligand)), units = c("cm"), dpi = 200 ) 
#part_name<-paste0(part_name,"din/")
if (!dir.exists(paste0("tost/"))) {dir.create(paste0("tost/"))}
if (!dir.exists(paste0("tcl/"))) {dir.create(paste0("tcl/"))}
if (!dir.exists(paste0("interaction_path/"))) {dir.create(paste0("interaction_path/"))}
i<-1
write.csv(df_structure,"path_structures.csv",row.names = F)
for (i in 1:nrow(df_structure)) {
   if (!dir.exists(paste0("tost/",df_structure$receptor[i],"_",df_structure$ligand[i]))) {
      dir.create(paste0("tost/",df_structure$receptor[i],"_",df_structure$ligand[i]))}
   receptor_name<-paste0(part_name,"receptor_start/",df_structure$receptor[i],".pdb")
   ligand_name<-paste0("str_fin/",df_structure$name.x[i])
   pdb_receptor<-read.pdb(receptor_name)
   pdb_ligand<-read.pdb(ligand_name)
   pdb_complex<-cat.pdb(pdb_receptor, pdb_ligand, rechain=TRUE)
   
   #  pdb<-read.pdb(paste0(part_name,"str_fin/",df_merge$name.x[i]))
   write.pdb(pdb_complex,paste0("tost/",df_structure$receptor[i],"_",df_structure$ligand[i],"/",df_structure$number[i],"_",df_structure$name.x[i]))
}
df_select<-df_structure%>%select(receptor,ligand)
df_select<-unique(df_select)
i<-1
for (i in 1:nrow(df_select)) {
   
   df_temp<-df_structure[df_structure$ligand==df_select$ligand[i]&df_structure$receptor==df_select$receptor[i],]
   df_temp<-df_temp%>%arrange(number)
   
   df_interactions<-read.csv(paste0( "interaction_surf/",df_temp$name.x[1],".csv"))
   df_interactions<-df_interactions%>%filter(persent_interactions==100)
   df_interactions<-df_interactions%>%mutate(name.x=df_temp$name.x[1])
   if (nrow(df_temp)>1){
      for (j in 2:nrow(df_temp)) {
         df_interactions_add<-read.csv(paste0( "interaction_surf/",df_temp$name.x[j],".csv"))
         df_interactions_add<-df_interactions_add%>%filter(persent_interactions==100)
         df_interactions_add<-df_interactions_add%>%mutate(name.x=df_temp$name.x[j])
         df_interactions<-rbind(df_interactions,df_interactions_add)
      }
   }
   
   df_interactions<-df_interactions%>%select(resno)
   df_interactions<-unique(df_interactions)
   df_temp<-df_temp%>%arrange((number))
   df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
   df_tcl[1,1]<-paste0('cd ', part_name,"din/tost/\n\n",
                       'color Display Background white\n',
                       'color Labels Atoms black\n',
                       'color Labels Bonds black\n\n',
                       'mol new {',paste0(part_name,"din/tost/",
                                          df_temp$receptor[1],"_",
                                          df_temp$ligand[1],"/",
                                          df_temp$number[1],"_",
                                          df_temp$name.x[1]),'} type {pdb}')
   for (j in 2:nrow(df_temp)) {
      df_tcl[j,1]<-paste0('mol addfile  {',paste0(part_name,"din/tost/",df_temp$receptor[j],"_",
                                                  df_temp$ligand[j],"/",
                                                  df_temp$number[j],"_",
                                                  df_temp$name.x[j]),'} type {pdb}')
      
   }
   
   df_tcl[j+1,1]<-paste0('mol modselect 0 ',0,' chain A\n',
                         'mol modmaterial 0 ',0,' Transparent\n',
                         'mol modstyle 0 ' ,0, ' NewCartoon\n',
                         'mol modselect 1 ',0,' chain B\n',
                         'mol addrep ',(0),'\n',
                         'mol modselect 1 ',0,' chain B\n',
                         'mol modstyle 1 ' ,0, ' Licorice\n')
   df_tcl[j+2,1]<-paste0('mol addrep ',(0),'\n',
                         'mol modselect 2 ',0,' chain A and resid ',paste(df_interactions$resno,collapse = " "),'\n',
                         'mol modstyle 2 ' ,0, ' NewCartoon\n')#,
   
   write.table(df_tcl,paste0("tcl/",df_select$receptor[i],"_",df_select$ligand[i],".tcl"),row.names = F,col.names = F,quote = F,sep = "\n",na="")
   df_interactions<-left_join(df_interactions,df_seq,by="resno")
   write.csv(df_interactions,file = paste0("interaction_path/",df_select$receptor[i],"_",df_select$ligand[i],".csv"),row.names = F)
   
}
