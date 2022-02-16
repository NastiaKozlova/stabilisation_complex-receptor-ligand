#!/usr/bin/env R
part_name = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)

setwd(part_name)
name<-list.files('start/structure/')

a<-c()
for(i in 1:length(name)){
  b<-strsplit(name[i],split = '.',fixed = T)[[1]][1]
  a<-c(a,b)
}
based_name<-a
name<-1
num_model<-100

for (name in 1:length(based_name)) {
  part<-paste0(part_name,based_name[name],'/MD/stabilisation/')

  setwd(paste0(part))
  if (!dir.exists(paste0('din/hbonds'))){dir.create(paste0('din/hbonds'))}
  if (!dir.exists(paste0('din/hbonds_log'))){dir.create(paste0('din/hbonds_log'))}
    df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
    df_tcl[1,1]<-paste('cd', part,"\npackage require hbonds")
    number_frame<-length(list.files(paste0('din/pdb_second/hbonds_',based_name[name])))
    if (number_frame>0){
      for (i in 0:(number_frame-1)) {
        df_tcl[i+2,1]<-paste0('mol new {protein/ionized_',based_name[name],'.psf} type {psf}')
        df_tcl[i+2,2]<-paste0('mol addfile {din/pdb_second/hbonds_',based_name[name],'/frame_',i,'.pdb} type {pdb}')
        df_tcl[i+2,3]<-paste0('set protein [atomselect top "protein" ]')
        df_tcl[i+2,4]<-paste0('set water [atomselect top "water" ]')
        df_tcl[i+2,5]<-paste0('hbonds -sel1 $protein -sel2 $water -writefile yes -upsel yes -frames all -dist 3.0 -ang 20 -plot no -outdir din -log hbonds_log/frame_',i,'.txt -writefile yes -outfile outfile -polar no -DA both -type all -detailout hbonds/frame_',i,'.txt')
        df_tcl[i+2,6]<-'mol delete all'
      }
      df_tcl[i+3,6]<-"\n\nexit now"
      write.table(df_tcl,file =paste0(part,'din/tcl/',based_name[name],'_hbond.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
      system(command = paste0('vmd -dispdev text -e ',part,'din/tcl/',based_name[name],'_hbond.tcl'),ignore.stdout=T,wait = T) 
    }
}
