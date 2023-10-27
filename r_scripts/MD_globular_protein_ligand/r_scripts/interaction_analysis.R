part_name = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
library(ggplot2)
library(cowplot)
#part_name<-part_start
setwd(part_name)
setwd("prepared_structures/")
name<-list.files("complex_structure/")
df_name<-data.frame(matrix(ncol=2,nrow = length(name)))
colnames(df_name)<-c("interactions_name","pdb_name")
df_name$interactions_name<-name
#a<-c()
i<-1
df_name<-df_name%>%mutate(ligand=NA)
for(i in 1:nrow(df_name)){
  b<-strsplit(df_name$interactions_name[i],split = ".",fixed = T)[[1]][1]
  b<-strsplit(b,split = "_",fixed = T)[[1]]
  df_name$pdb_name[i]<-paste0(b[1:length(b)],collapse = "_")
  df_name$ligand[i]<-paste0(b[2],collapse = "_")
  #  a<-c(a,b)
}
df_name<-df_name%>%select(pdb_name,ligand)
df_name<-df_name%>%filter(ligand!="DMPE")
df_name<-unique(df_name)
#if(!dir.exists("complex_structure_center")){dir.create("complex_structure_center")}

if(!dir.exists("interactions_center")){dir.create("interactions_center")}
if(!dir.exists("atom_interactions_center")){dir.create("atom_interactions_center")}
if(!dir.exists("make_picture_tcl_center")){dir.create("make_picture_tcl_center")}
#df_merge<-df_merge%>%mutate(complex_name=paste0(receptor,"_",ligand,"_",size_of_group))
#i<-1
#df_merge<-df_merge%>%mutate(tested=F)
j<-1
for (j in 1:nrow(df_name)) {
  if(!file.exists(paste0("interactions_center/",df_name$pdb_name[j],".csv"))){
    #      file.exists(paste0("complex_structure/",df_name$pdb_name[j],".pdb.pdb"))
    pdb<-read.pdb(paste0("complex_structure/",df_name$pdb_name[j],".pdb.pdb"))
    protein.int<-atom.select(pdb,"protein")
    ligand.int<-atom.select(pdb, "protein", inverse=TRUE)
    protein <- trim.pdb(pdb, protein.int)
    ligand <- trim.pdb(pdb, ligand.int)
    #      b<-read.pdb(paste0("pdb_second/",df_name$pdb_name[j],"/frame_",df_all$new_number[j],".pdb"))
    bs<-binding.site(protein,ligand)
    m<-bs$resnames
    a<-c()
    b<-c()
    y<-1
    for (y in 1:length(m)) {
      p<-strsplit(m[y],split = " ",fixed = T)[[1]][2]
      a<-c(a,p)
      p<-strsplit(m[y],split = " ",fixed = T)[[1]][1]
      b<-c(b,p)
    }
    a<-as.numeric(a)
    df_protein<-data.frame(matrix(ncol=2,nrow=length(a)))
    colnames(df_protein)<-c("resno","resid")
    df_protein$resid<-b
    df_protein$resno<-a
    write.csv(df_protein,
              paste0("interactions_center/",df_name$pdb_name[j],".csv"),
              row.names = F)
  }
}
j<-1
df_protein<-read.csv(paste0("interactions_center/",df_name$pdb_name[1],".csv"),
                     stringsAsFactors = F)
df_protein<-df_protein%>%mutate(pdb_name=df_name$pdb_name[1])
if(nrow(df_name)>1){
  for (j in 2:nrow(df_name)) {
    df_protein_add<-read.csv(paste0("interactions_center/",df_name$pdb_name[j],".csv"),
                             stringsAsFactors = F)
    df_protein_add<-df_protein_add%>%mutate(pdb_name=df_name$pdb_name[j])
    df_protein<-rbind(df_protein,df_protein_add)
  }
}
write.csv(df_protein,file = paste0("interactions_center_combined.csv"),row.names = F)

for (j in 1:nrow(df_name)) {
  
  df_interactions<-read.csv(paste0("interactions_center/",df_name$pdb_name[j],".csv"),stringsAsFactors = F)
  
  df_interactions<-df_interactions%>%select(resid,resno)
  df_interactions<-unique(df_interactions)
  pdb<-read.pdb(paste0("complex_structure/",df_name$pdb_name[j],".pdb.pdb"))
  protein.int<-atom.select(pdb,"protein")
  ligand.int<-atom.select(pdb, "protein", inverse=TRUE)
  protein <- trim.pdb(pdb, protein.int)
  ligand <- trim.pdb(pdb, ligand.int)
  #   v_binding_test<-binding.site(protein,ligand,cutoff = 12)$resno
  df_protein<-protein$atom
  df_ligand<-ligand$atom
  df_protein<-df_protein[df_protein$resno%in%df_interactions$resno,]
  for (q in 1:nrow(df_protein)) {
    df_protein$alt[q]<-strsplit(df_protein$elety[q],split = "",fixed = T)[[1]][1]
  }
  for (q in 1:nrow(df_ligand)) {
    df_ligand$alt[q]<-strsplit(df_ligand$elety[q],split = "",fixed = T)[[1]][1]
  }
  df_test<-full_join(df_protein,df_ligand,by="type",
                     relationship = "many-to-many")
  df_test<-df_test%>%mutate(length=sqrt((x.x-x.y)^2+(y.x-y.y)^2+(z.x-z.y)^2))
  df_test<-df_test%>%filter(length<12)
  df_test<-df_test%>%filter(alt.x!=alt.y)
  df_test<-df_test%>%select(eleno.x,  elety.x,  alt.x,    resid.x,  resno.x,x.x,y.x,z.x,
                            eleno.y,  elety.y,  alt.y,    resid.y,  resno.y,x.y,y.y,z.y, length)
  df_interaction<-df_test%>%filter(elety.x!="CA")
  df_interaction<-df_interaction%>%filter(elety.x!="N")
  df_interaction<-df_interaction%>%group_by(eleno.y)%>%mutate(length_test=min(length))
  df_interaction<-df_interaction%>%group_by(eleno.y)%>%filter(length_test==length)
  df_interaction<-ungroup(df_interaction)
  df_interaction<-df_interaction%>%group_by(resno.x)%>%mutate(length_test=min(length))
  df_interaction<-df_interaction%>%group_by(resno.x)%>%filter(length_test==length)
  df_interaction<-ungroup(df_interaction)
  write.csv(df_interaction,file = paste0("atom_interactions_center/",df_name$pdb_name[j],".csv"),row.names = F)
}

#df_name<-df_name%>%filter(ligand!="DMPE")
for (j in 1:nrow(df_name)) {
  
  df_interactions<-read.csv(paste0("atom_interactions_center/",df_name$pdb_name[j],".csv"),stringsAsFactors = F)
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  df_tcl[1,1]<-paste0('cd ', part_name,"prepared_structures/complex_structure/\n\n",
                      'mol new {',df_name$pdb_name[j],'.pdb.psf} type {psf}\n',
                      'mol addfile {',df_name$pdb_name[j],'.pdb.pdb} type {pdb}\n')
  df_tcl[1,2]<-paste0('color Display Background white\n',
                      'color Labels Atoms black\n',
                      'color Labels Bonds black\n\n')
  df_tcl[1,3]<-paste0('mol modselect 0 ',j-1,' protein\n',
                      'mol modmaterial 0 ',(j-1),' Transparent\n',
                      'mol modstyle 0 ' ,j-1, ' NewCartoon\n')#,
  df_tcl[1,4]<-paste0('mol addrep ',(j-1),'\n',
                      'mol modselect 1 ',j-1,' not protein\n',
                      'mol modmaterial 1 ',(j-1),' Opaque\n',
                      'mol modstyle 1 ' ,j-1, ' Licorice\n')#,
  b<-paste('resid ',paste0(unique(df_interactions$resno.x),collapse = " "))
  df_tcl[1,5]<-paste0('mol addrep ',(j-1),'\n',
                      'mol modselect 2 ',j-1,' protein and ',b,'\n',
                      'mol modmaterial 2 ',(j-1),' Opaque\n',
                      'mol modstyle 2 ' ,j-1, ' Licorice\n',
                      'mol modcolor 2 ',(j-1),' ResType\n')#,
  
  b<-paste('(resid ',df_interaction$resno.x,' and name ',df_interaction$elety.x," and resname ",df_interaction$resid.x,")")
  a<-paste0(b,collapse = " or ")
  df_tcl[1,6]<-paste0('set all [atomselect ',(j-1),' "',a,'"]\n',
                      'set i ',(j-1),"\n",
                      'foreach atom [$all list] {\n',
                      '  label add Atoms ',(j-1),'/$atom\n',
                      '  incr i\n}\n',
                      '$all delete\n\n')
  
  for (p in 1:nrow(df_interactions)) {
    df_tcl[(p+1),1]<-paste0('set atomID1 [[atomselect ',(j-1),' "(resid ',df_interactions$resno.x[p],
                            ' and name ',df_interactions$elety.x[p],')"] list]\n')
    df_tcl[(p+1),2]<-paste0('set atomID2 [[atomselect ',(j-1),
                            ' "(x > ',df_interactions$x.y[p]-0.5,' and x < ',df_interactions$x.y[p]+0.5,
                            ' and y > ',df_interactions$y.y[p]-0.5,' and y < ',df_interactions$y.y[p]+0.5,
                            ' and z > ',df_interactions$z.y[p]-0.5,' and z < ',df_interactions$z.y[p]+0.5,')"] list]\n')
    
    df_tcl[(p+1),3]<-paste0('label add Bonds ',(j-1),'/$atomID1 ',(j-1),'/$atomID2\n')
  }
  df_tcl[is.na(df_tcl)]<-""
  write.csv(df_tcl,paste0("make_picture_tcl_center/",df_name$pdb_name[j],".tcl"),row.names = F)
}

#df_merge<-df_merge%>%filter(tested)
df_tcl<-read.csv(paste0("make_picture_tcl_center/",df_name$pdb_name[1],".tcl"),stringsAsFactors = F)
i<-2
if(nrow(df_name)>1){
  for (i in 2:nrow(df_name)) {
    df_tcl_add<-read.csv(paste0("make_picture_tcl_center/",df_name$pdb_name[i],".tcl"),stringsAsFactors = F)
    df_tcl<-rbind(df_tcl,df_tcl_add)
  }
}
write.table(df_tcl,paste0("make_picture_tcl_center.tcl"),row.names = F,col.names = F,quote = F,sep = "\n",na="")
