part_name = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
library(readr)

part_name<-paste0(part_name,"docking_first")
setwd(part_name)
name<-"start"

part<-1
sort_hbonds<-function(file_name,frame_number){
  df_hbonds<-read_tsv(file_name,skip = 2,col_names = F)
  colnames(df_hbonds)<-c("donor","acceptor","occupancy")
  df_hbonds<-df_hbonds%>%mutate(frame=frame_number)
  df_hbonds<-df_hbonds%>%mutate(amino=NA)
  df_hbonds<-df_hbonds%>%mutate(number=NA)
  df_hbonds<-df_hbonds%>%mutate(type=NA)
  for (i in 1:nrow(df_hbonds)) {
    dd<-strsplit(df_hbonds$donor[i],split = "-",fixed = T)[[1]][1]
    dd<-strsplit(dd,split = "",fixed = T)[[1]]
    number_d<-paste0(dd[4:length(dd)],collapse = "")
    amino_d<-paste0(dd[1:3],collapse = "")
    aa<-strsplit(df_hbonds$acceptor[i],split = "-",fixed = T)[[1]][1]
    aa<-strsplit(aa,split = "",fixed = T)[[1]]
    number_a<-paste0(aa[4:length(aa)],collapse = "")
    amino_a<-paste0(aa[1:3],collapse = "")
    if (amino_a=="wat") {
      df_hbonds$type[i]<-"donor"
      df_hbonds$number[i]<-number_d
      df_hbonds$amino[i]<-amino_d
    }
    if (amino_d=="wat") {
      df_hbonds$type[i]<-"acceptor"
      df_hbonds$number[i]<-number_a
      df_hbonds$amino[i]<-amino_a
    }
  }
  df_hbonds<-df_hbonds%>%select(number,amino,number,type,frame)
  df_hbonds<-unique(df_hbonds)
  df_hbonds<-df_hbonds%>%mutate(number=as.numeric(number))
  return(df_hbonds)
}

if (!dir.exists(paste0('hbonds_test'))) {dir.create(paste0('hbonds_test'))}

  #prepare psf and pdb to run NAMD
  
  df_psfgen<-data.frame(matrix(ncol = 1,nrow = 1))
  df_psfgen[1,1]<-paste0('cd ',part_name,'\n',
  'mol delete all\n',
  'package require psfgen \n',
  'resetpsf\n',
  'topology toppar/top_all36_prot.rtf\n',
  'topology toppar/toppar_water_ions_namd.str\n',
  
  'pdbalias residue HIS HSE\n',
  'pdbalias atom ILE CD1 CD\n',
  'segment U { pdb receptor_start/start.pdb\n',
  '}\n',
  '\ncoordpdb receptor_start/start.pdb U\n',
  'regenerate angles dihedrals\n',
  'guesscoord\n',
  'writepdb hbonds_test/',name,'.pdb\n',
  'writepsf hbonds_test/',name,'.psf\n',
  'mol delete all\n\nexit now')
  write.table(df_psfgen,paste0('hbonds_test/psfgen_',name,'.tcl'),col.names = F,row.names = F,quote = F)
  system(command = paste0("vmd -dispdev text -e ",part_name,'/hbonds_test/psfgen_',name,'.tcl'),ignore.stdout=T,wait = T) 
  print(paste0(name))
  pdb_start<-read.pdb(paste0('hbonds_test/',name,'.pdb'))
  df_pdb<-pdb_start$atom
  
  x_min<-min(df_pdb$x)-20
  y_min<-min(df_pdb$y)-20
  z_min<-min(df_pdb$z)-20
  
  x_max<-max(df_pdb$x)+20
  y_max<-max(df_pdb$y)+20
  z_max<-max(df_pdb$z)+20
  #prepare script to solvate and ionize protein pstructure
  df_psfgen<-data.frame(matrix(ncol = 1,nrow = 1))
  df_psfgen[1,1]<-paste0('cd ',part_name,'/hbonds_test\n',
                         'package require solvate \n','package require autoionize \n',
                         'solvate ',name,'.psf ',name,'.pdb -o solvate_',name,' -b 1.5 -minmax {{',x_min,' ',y_min,' ',z_min,'} {',x_max,' ',y_max,' ',z_max,'}}\n',
                         'autoionize -psf solvate_',name,'.psf -pdb solvate_',name,'.pdb -sc 0.15 -o ionized_',name,
                         '\nmol delete all\n\nexit now')
  write.table(df_psfgen,paste0('hbonds_test/solvate_',name[part],'.tcl'),col.names = F,row.names = F,quote = F)
  system(command = paste0("vmd -dispdev text -e ",part_name,'/hbonds_test/solvate_',name[part],'.tcl'),ignore.stdout=T,wait = T) 
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  df_tcl[1,1]<-paste('cd', part_name,"\npackage require hbonds")
  df_tcl[1,2]<-paste0('mol new {hbonds_test/ionized_',name,'.psf} type {psf}')
  df_tcl[1,3]<-paste0('mol addfile {hbonds_test/ionized_',name,'.pdb} type {pdb}')
  df_tcl[1,4]<-paste0('set protein [atomselect top "protein" ]')
  df_tcl[1,5]<-paste0('set water [atomselect top "water" ]')
  df_tcl[1,6]<-paste0('hbonds -sel1 $protein -sel2 $water -writefile yes -upsel yes -frames all -dist 3.0 -ang 20 -plot no -outdir hbonds_test -log hbonds_log.txt -writefile yes -outfile outfile -polar no -DA both -type all -detailout hbonds_test.txt')
  df_tcl[1,7]<-'mol delete all'
  df_tcl[1,8]<-"\n\nexit now"
  
  write.table(df_tcl,file =paste0(part_name,'/hbonds_test/',name,'_hbond.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)

  system(command = paste0("vmd -dispdev text -e ",part_name,'/hbonds_test/',name,'_hbond.tcl'),ignore.stdout=T,wait = T) 
  file_name<-paste0("hbonds_test/hbonds_test.txt")
  df_bonds<-sort_hbonds(file_name=file_name,frame_number=1)
  pdb<-read.pdb(paste0("receptor_start/",name,".pdb"))
  df_pdb<-pdb$atom
  df_pdb<-df_pdb%>%filter(elety=="CA")
  df_pdb<-df_pdb%>%select(resid,resno)
  df_bonds<-left_join(df_bonds,df_pdb,c("number"="resno","amino"="resid"))
  df_bonds<-unique(df_bonds)
  df_bonds<-df_bonds%>%mutate(persent=100)
  write.csv(df_bonds,paste0("hbonds.csv"),row.names = F)  
  