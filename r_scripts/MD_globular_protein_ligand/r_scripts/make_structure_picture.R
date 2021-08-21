part_start <- commandArgs(trailingOnly=TRUE)
#part_start<-part_name
setwd(part_start)

library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)

setwd(part_start)
if(file.exists(paste0("complex.csv"))){
  df_complex<-read.csv(paste0("complex.csv"),stringsAsFactors = F)
  v_ligand<-unique(df_complex$ligand)
  i<-1
  j<-3
  for (i in 1:length(v_ligand)) {
    df_complex_TEMP<-df_complex%>%filter(ligand==v_ligand[i])
    a<-0
    for (j in 1:nrow(df_complex_TEMP)) {
      if(file.exists(paste0(df_complex_TEMP$complex_name[j],"/MD/stabilisation/quench/quench_",df_complex_TEMP$complex_name[j],".dcd")))
        a<-a+1  
    }
    
    if(a==nrow(df_complex_TEMP)){
      for (j in 1:nrow(df_complex_TEMP)) {
        
      
      df_interaction<-read.csv(paste0('fin_structure/interactions_atom/',df_complex_TEMP$complex_name[j],".csv"),stringsAsFactors = F)
      df_interaction<-df_interaction%>%filter(Total<(-1))
      
      df_interaction<-df_interaction%>%mutate(sort=paste(resno.x,elety.y))
      df_interaction<-df_interaction%>%group_by(sort)%>%mutate(min_Total=min(Total))
      df_interaction<-df_interaction%>%filter(min_Total==Total)
      
      df_interaction<-df_interaction%>%group_by(resno.x)%>%mutate(min_Total=min(Total))
      df_interaction<-df_interaction%>%filter(min_Total==Total)
      
      df_interactions<-df_interaction%>%select(resid.x,resno.x,elety.x,resid.y,resno.y,elety.y)

      df_interactions_test<-df_interaction#<-rbind(df_interactions,df_interactions_add)
                                                    #length,elec,vdw,nonb,Total,bond_type)
      df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
      df_tcl[1,1]<-paste0('cd ', part_start,"\n")
      df_tcl[1,2]<-paste0('mol new {fin_structure/',v_ligand[i],"/",df_complex_TEMP$complex_name[j],'.psf} type {psf}')
      df_tcl[1,3]<-paste0('mol addfile {fin_structure/',v_ligand[i],"/",df_complex_TEMP$complex_name[j],'.pdb} type {pdb}')
      b<-paste('(resid ',df_interactions_test$resno.x,' and name ',df_interactions_test$elety.x," and resname ",df_interactions_test$resid.x,")")
      a<-paste0(b,collapse = " or ")
      df_tcl[1,4]<-paste0('set all [atomselect ',(j-1),' "',a,'"]')
      df_tcl[1,5]<-paste0('set i ',(j-1))
      df_tcl[1,6]<-paste0('foreach atom [$all list] {')
      df_tcl[1,7]<-paste0('  label add Atoms ',(j-1),'/$atom')
      df_tcl[1,8]<-paste0('  incr i\n}')
      df_tcl[1,9]<-paste0('$all delete')
      for (p in 1:nrow(df_interaction)) {
        
      
        df_tcl[(p+1),1]<-paste0('set atomID1 [[atomselect ',(j-1),' "(resid ',df_interaction$resno.x[p],
                            ' and name ',df_interaction$elety.x[p],
                            " and resname ",df_interaction$resid.x[p],')"] list]')
        df_tcl[(p+1),2]<-paste0('set atomID2 [[atomselect ',(j-1),' "(resid ',df_interaction$resno.y[p],
                            ' and name ',df_interaction$elety.y[p],
                            " and resname ",df_interaction$resid.y[p],')"] list]')
        
        df_tcl[(p+1),3]<-paste0('label add Bonds ',(j-1),'/$atomID1 ',(j-1),'/$atomID2')
      }
      df_tcl[is.na(df_tcl)]<-""
      write.table(df_tcl,paste0("fin_structure/make_picture_tcl/",df_complex_TEMP$complex_name[j],"_test.tcl"),row.names = F,col.names = F,quote = F,sep = "\n",na="")
      print(paste0(df_complex_TEMP$structure_name[j],' resid ',paste0(unique(df_interaction$resno.x),collapse = " ")))
      }
    }
  }
}
      
