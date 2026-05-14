part_analysis <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(dplyr)
library(ggplot2)
#paste0(part_analysis)
df_all<-read.csv(paste0(part_analysis,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%select(receptor,ligand)
df_all<-unique(df_all)
part_name<-paste0(part_analysis,"din/")
setwd(part_name)
i<-1
#for (i in 1:nrow(df_all)) {
#  if(!file.exists(paste0("interaction_fin/",df_all$receptor[i],"_",df_all$ligand[i],".csv"))){
#    df_all$receptor[i]<-NA
#  }
#}
#df_all<-df_all%>%filter(!is.na(receptor))
df_merge<-read.csv(paste0(part_name,"log_merge.csv"),stringsAsFactors = F)
#df_merge<-semi_join(df_merge,df_all)
df_merge<-df_merge%>%select(name.x,receptor,ligand)
df_merge<-unique(df_merge)
i<-1
#if (dir.exists(paste0("complex_structure/"))) {system(command = paste0("rm -r ",part_analysis,"din/complex_structure/"))}
if (dir.exists(paste0("make_picture_tcl/"))) {system(command = paste0("rm -r ",part_analysis,"din/make_picture_tcl/"))}
#if(!dir.exists("complex_structure")){dir.create("complex_structure")}
if(!dir.exists("make_picture_tcl")){dir.create("make_picture_tcl")}
df_merge<-df_merge%>%mutate(complex_name=paste0(receptor,"_",ligand))
i<-1
for (i in 1:nrow(df_merge)) {
  
  if(!file.exists(paste0("interaction_merged/",df_merge$name.x[i]))){
    df_merge$name.x[i]<-NA
  }
}
df_merge<-df_merge%>%filter(!is.na(name.x))
df_merge<-df_merge%>%mutate(number=1:nrow(df_merge)-1)
for (i in 1:nrow(df_merge)) {
  #  df_hbonds<-read.csv(paste0(part_analysis,"/hbonds.csv"),stringsAsFactors = F)
  #  df_hbonds<-df_hbonds%>%filter(persent>0)
  #  if(file.exists(paste0("interaction_merged/",df_merge$name.x[i]))){
  df_interactions<-read.csv(paste0("interaction_merged/",df_merge$name.x[i]),stringsAsFactors = F)
  #df_interactions<-df_interactions%>%filter(persent_interactions>50)
  df_interactions<-df_interactions%>%filter(!is.na(Motif))
  df_interactions<-df_interactions%>%select(resid,resno,Motif)
  df_interactions<-unique(df_interactions)
  
  receptor_name<-paste0(part_analysis,"receptor_start/",df_merge$receptor[i],".pdb")
  ligand_name<-paste0(part_name,"str_fin/",df_merge$name.x[i])
  protein<-read.pdb(receptor_name)
  ligand<-read.pdb(ligand_name)
  df_ligand<-ligand$atom
  for (q in 1:nrow(df_ligand)) {
    df_ligand$alt[q]<-strsplit(df_ligand$elety[q],split = "",fixed = T)[[1]][1]
  }
  df_ligand<-df_ligand%>%filter(alt!="C")
  df_protein<-protein$atom
  df_protein<-df_protein[df_protein$resno%in%df_interactions$resno,]
  for (q in 1:nrow(df_protein)) {
    df_protein$alt[q]<-strsplit(df_protein$elety[q],split = "",fixed = T)[[1]][1]
  }
  df_protein<-df_protein%>%filter(alt!="C")
  df_protein<-df_protein%>%filter(alt!="H")
  df_protein<-df_protein%>%filter(elety!="N")
  df_protein<-df_protein%>%filter(elety!="O")
  df_protein_start<-left_join(df_protein,df_ligand,by="type",
                              relationship = "many-to-many")
#  df_protein_start<-df_protein_start%>%filter(resno.x<0)
#  for (l in 1:nrow(df_ligand)) {
#    ligand.int<-atom.select(ligand,eleno=df_ligand$eleno[l])
#    ligand.pdb <- trim.pdb(ligand, ligand.int)
#    df_ligand.pdb<-df_ligand[l,]
    #v_binding_test<-binding.site(protein,ligand.pdb,cutoff=5)$resno
    #if(length(v_binding_test)>0){
    #df_protein_test<-df_protein[df_protein$resno%in%v_binding_test,]
#    df_protein_test<-left_join(df_protein,df_ligand.pdb,by="type",
#                               relationship = "many-to-many")
#    df_protein_start<-rbind(df_protein_start,df_protein_test)
#    print(paste(l,nrow(df_protein_test)))
    #}
#  }
  df_protein_start<-unique(df_protein_start)

  df_test<-df_protein_start%>%mutate(length=sqrt((x.x-x.y)^2+(y.x-y.y)^2+(z.x-z.y)^2))
  #df_test<-df_test%>%filter(length<12)
  
  df_test<-df_test%>%filter(alt.x!=alt.y)
  
  
  
  df_test<-df_test%>%select(eleno.x,  elety.x,  alt.x,    resid.x,  resno.x,x.x,y.x,z.x,
                            eleno.y,  elety.y,  alt.y,    resid.y,  resno.y,x.y,y.y,z.y, length)
  df_interaction<-df_test%>%group_by(resno.x)%>%mutate(length_test=min(length))
  print(nrow(df_interaction))
  
  df_interaction<-df_interaction%>%filter(length_test==length) 
  print(nrow(df_interaction))
  
  df_interaction<-df_interaction%>%filter(alt.y!="C")
  df_interaction<-df_interaction%>%filter(alt.x!="C")
  df_interaction<-df_interaction%>%filter(alt.x!="H")
  #df_interaction<-df_interaction%>%filter(elety.x!="N")
  #df_interaction<-df_interaction%>%filter(elety.x!="O")
  
  df_interaction<-ungroup(df_interaction)
  
  df_interaction<-df_interaction%>%group_by(eleno.x)%>%
    mutate(length_test=min(length))%>%
    filter(length_test==length)
  df_interaction<-df_interaction%>%group_by(eleno.y)%>%filter(length_test<12)
  #df_interaction<-ungroup(df_interaction)
  
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  df_tcl[1,1]<-paste0('cd ', part_name,"complex_structure/\n\n",
                      'mol new {',df_merge$name.x[i],'} type {pdb}')
  b<-paste('(resid ',df_interaction$resno.x,' and name ',df_interaction$elety.x," and resname ",df_interaction$resid.x,")")
  a<-paste0(b,collapse = " or ")
  df_tcl[1,2]<-paste0('set all [atomselect ',(df_merge$number[i]),' "',a,'"]\n',
                      'set i ',(df_merge$number[i]),"\n",
                      'foreach atom [$all list] {\n',
                      '  label add Atoms ',(df_merge$number[i]),'/$atom\n',
                      '  incr i\n}\n',
                      '$all delete\n\n',
                      'color Display Background white\n',
                      'color Labels Atoms black\n',
                      'label textsize 1\n',
                      'label textthickness 3\n',
                      'color Labels Bonds black\n\n')
  
  df_tcl[1,3]<-paste0('mol modselect 0 ',df_merge$number[i],' protein\n',
                      'mol modmaterial 0 ',(df_merge$number[i]),' Transparent\n',
                      'mol modstyle 0 ' ,df_merge$number[i], ' NewCartoon\n')
  
  if (nrow(df_interaction)>0){
    df_tcl[1,4]<-paste0('mol selection resname ',paste0(unique(df_interaction$resid.y),collapse = " "))
    df_tcl[1,5]<-paste0('mol modmaterial 1 ',(df_merge$number[i]),' Opaque\n',
                        'mol addrep ',(df_merge$number[i]),'\n',
                        'mol modstyle 1 ',(df_merge$number[i]),' Licorice \n',
                        'mol modcolor 1 ',(df_merge$number[i]),' Name\n')
    df_tcl[1,6]<-paste0('mol selection (resid ',paste0(unique(df_interactions$resno),collapse = " "),")")
    df_tcl[1,7]<-paste0(' mol modmaterial 2 ',(df_merge$number[i]),' Opaque\n',
                        'mol addrep ',(df_merge$number[i]),'\n',
                        'mol modstyle 2 ',(df_merge$number[i]),' Licorice 0.1 12 12 \n',
                        'mol modcolor 2 ',(df_merge$number[i]),' ColorID 0')
    df_tcl[1,8]<-paste0('mol selection (resid ',paste0(unique(df_interaction$resno.x),collapse = " "),")")
    df_tcl[1,9]<-paste0(' mol modmaterial 3 ',(df_merge$number[i]),' Opaque\n',
                        'mol addrep ',(df_merge$number[i]),'\n',
                        'mol modstyle 3 ',(df_merge$number[i]),' Licorice 0.1 12 12 \n',
                        'mol modcolor 3 ',(df_merge$number[i]),' ColorID 1')
    df_interactions_motif<-df_interactions%>%filter(!is.na(Motif))
    v_center<-df_interactions_motif$resno[df_interactions_motif$resno%in%df_interaction$resno.x]
    v_center<-paste0(v_center,collapse = " ")
    df_tcl[1,10]<-paste0('mol selection (resid ',v_center,")")
    df_tcl[1,11]<-paste0(' mol modmaterial 4 ',(df_merge$number[i]),' Opaque\n',
                        'mol addrep ',(df_merge$number[i]),'\n',
                        'mol modstyle 4 ',(df_merge$number[i]),' Licorice 0.1 12 12 \n',
                        'mol modcolor 4 ',(df_merge$number[i]),' Type')

    for (p in 1:nrow(df_interaction)) {
      df_tcl[(p+1),1]<-paste0('set sel1 [atomselect ',(df_merge$number[i]),' "(resid ',df_interaction$resno.x[p],
                              ' and name ',df_interaction$elety.x[p],
                              " and resname ",df_interaction$resid.x[p],')"]')
      df_tcl[(p+1),2]<-paste0('set sel2 [atomselect ',(df_merge$number[i]),' "(resid ',df_interaction$resno.y[p],
                              ' and name ',df_interaction$elety.y[p],
                              " and resname ",df_interaction$resid.y[p],')"]')
      #        df_tcl[(p+1),2]<-paste0('set sel2 [atomselect ',(df_merge$number[i]),
      #                                ' "(x > ',df_interaction$x.y[p]-0.5,' and x < ',df_interaction$x.y[p]+0.5,
      #                                ' and y > ',df_interaction$y.y[p]-0.5,' and y < ',df_interaction$y.y[p]+0.5,
      #                                ' and z > ',df_interaction$z.y[p]-0.5,' and z < ',df_interaction$z.y[p]+0.5,')"]')
      
      df_tcl[(p+1),3]<-paste0('lassign [$sel1 get {x y z}] pos1\n',
                              'lassign [$sel2 get {x y z}] pos2\n',
                              'draw color black\n',
                              'draw line $pos1 $pos2 width 3')
      #                                'draw line $pos1 $pos2 style dotted width 2')
      
      
      # make the two selections
      #        set sel1 [atomselect top "resid 45 and name CA"]
      #        set sel2 [atomselect top "resid 99 and name N"]
      # get the coordinates
      #        lassign [$sel1 get {x y z}] pos1
      #        lassign [$sel2 get {x y z}] pos2
      # draw a white line between the two atoms
      #        draw color white
      #        draw line $pos1 $pos2 style dotted width 2
      
    }
  }
  df_tcl[is.na(df_tcl)]<-""
  write.csv(df_tcl,paste0("make_picture_tcl/",df_merge$name.x[i],".tcl"),row.names = F)
  #  }else{print(i)}
} 

#v_structure<-list.files("make_picture_tcl/")
df_tcl<-read.csv(paste0("make_picture_tcl/",df_merge$name.x[1],".tcl"),stringsAsFactors = F)
write.table(df_tcl,paste0("make_picture_tcl.tcl"),row.names = F,col.names = F,quote = F,sep = "\n",na="")
i<-2
for (i in 2:nrow(df_merge)) {
  df_tcl_add<-read.csv(paste0("make_picture_tcl/",df_merge$name.x[i],".tcl"),stringsAsFactors = F)
  df_tcl<-rbind(df_tcl,df_tcl_add)
}
#write.table(df_tcl,paste0("make_picture_tcl.tcl"),row.names = F,col.names = F,quote = F,sep = "\n",na="")
