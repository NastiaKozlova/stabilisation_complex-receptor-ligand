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
  j<-1
  for (i in 1:length(v_ligand)) {
    df_complex_TEMP<-df_complex%>%filter(ligand==v_ligand[i])
    a<-0
    for (j in 1:nrow(df_complex_TEMP)) {
      if(file.exists(paste0(df_complex_TEMP$complex_name[j],"/MD/stabilisation/quench/quench_",df_complex_TEMP$complex_name[j],".dcd")))
        a<-a+1  
    }
    
    if(a==nrow(df_complex_TEMP)){
      df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
      df_tcl[1,1]<-paste0('cd ', part_start,'\n')
      
      if(!dir.exists("fin_structure")){dir.create("fin_structure")}
      if(!dir.exists("fin_structure/tcl")){dir.create("fin_structure/tcl")}
      if(!dir.exists("fin_structure/plot")){dir.create("fin_structure/plot")}
      if(!dir.exists("fin_structure/interactions")){dir.create("fin_structure/interactions")}
      if(!dir.exists("fin_structure/interactions_atom")){dir.create("fin_structure/interactions_atom")}
      
      if(!dir.exists(paste0("fin_structure/",v_ligand[i]))){dir.create((paste0("fin_structure/",v_ligand[i])))}
      if(!dir.exists(paste0("fin_structure/",v_ligand[i],"/full_structure"))){dir.create((paste0("fin_structure/",v_ligand[i],"/full_structure")))}
      for (j in 1:nrow(df_complex_TEMP)) {
        if(!dir.exists(paste0("fin_structure/",v_ligand[i],"/",df_complex_TEMP$complex_name[j]))){dir.create(paste0("fin_structure/",v_ligand[i],"/",df_complex_TEMP$complex_name[j]))}
        df_interactions<-read.csv(paste0(df_complex_TEMP$complex_name[j],"/MD/stabilisation/din/intercations.csv"),stringsAsFactors = F)
        df_interactions<-df_interactions%>%filter(number<1000)
        df_interactions<-df_interactions%>%filter(Total==min(df_interactions$Total))
        
        #Secondary structure counting
#        df_tcl<-data.frame(matrix(nrow = 1,ncol = 7))
        df_tcl[(j+1),1]<-paste0('mol new {',df_complex_TEMP$complex_name[j],'/MD/stabilisation/protein/ionized_',df_complex_TEMP$complex_name[j],'.psf} type {psf}')
        df_tcl[(j+1),2]<-paste0('mol addfile {',df_complex_TEMP$complex_name[j],'/MD/stabilisation/quench/quench_',df_complex_TEMP$complex_name[j],'.dcd} type {dcd} first ',df_interactions$number[1],' last ',df_interactions$number[1],' step 1 waitfor all')
        df_tcl[(j+1),3]<-paste0('set nf [molinfo top get numframes]')
        df_tcl[(j+1),4]<-paste0('[atomselect top all] writepdb fin_structure/',v_ligand[i],'/full_structure/',df_complex_TEMP$complex_name[j],'.pdb\n')
        df_tcl[(j+1),5]<-paste0('[atomselect top all] writepsf fin_structure/',v_ligand[i],'/full_structure/',df_complex_TEMP$complex_name[j],'.psf\n')
        df_tcl[(j+1),6]<-'mol delete all\n\n'
      }
      df_tcl[(j+2),1]<-paste0('mol delete all \n \n \n exit now')
      write.table(df_tcl,paste0('fin_structure/tcl/full_structure_',v_ligand[i],'.tcl'),na = "\n",col.names = F,row.names = F,quote = F,sep = "\n")
      system(command = paste0("vmd -dispdev text -e ",part_start,'fin_structure/tcl/full_structure_',v_ligand[i],'.tcl'),ignore.stdout=T,wait = T) 
      j<-1
      for (j in 1:nrow(df_complex_TEMP)) {
        df_interactions<-read.csv(paste0(df_complex_TEMP$complex_name[j],"/MD/stabilisation/din/intercations.csv"),stringsAsFactors = F)
        df_interactions<-df_interactions%>%filter(number<1000)
        df_interactions<-df_interactions%>%filter(Total==min(df_interactions$Total))
        df_interaction<-read.csv(paste0(df_complex_TEMP$complex_name[j],'/MD/stabilisation/din/interaction/',df_complex_TEMP$complex_name[j],'/frame_',df_interactions$number[1],'.csv'),stringsAsFactors = F)
#        print(paste0(df_complex_TEMP$center[j]," resid ",paste0(df_interaction$resid,collapse = " ")))
        pdb<-read.pdb(paste0(df_complex_TEMP$complex_name[j],'/MD/stabilisation/din/pdb_second_complex/',df_complex_TEMP$complex_name[j],'/frame_',df_interactions$number[1],'.pdb'))
        receptor.int<-atom.select(pdb,"protein")
        ligand.int<-atom.select(pdb,"ligand")
        protein<-trim(pdb,receptor.int)
        ligand<-trim(pdb,ligand.int)
        df_protein<-protein$atom
        df_ligand<-ligand$atom
        df_protein<-df_protein[df_protein$resno%in%df_interaction$resid,]
        q<-1
        for (q in 1:nrow(df_protein)) {
          df_protein$alt[q]<-strsplit(df_protein$elety[q],split = "",fixed = T)[[1]][1]
        }
        for (q in 1:nrow(df_ligand)) {
          df_ligand$alt[q]<-strsplit(df_ligand$elety[q],split = "",fixed = T)[[1]][1]
        }
        df_ligand<-df_ligand%>%filter(alt!="C")
        df_protein<-df_protein%>%filter(alt!="C")
        
        df_test<-full_join(df_protein,df_ligand,by="type")
        df_test<-df_test%>%mutate(length=sqrt((x.x-x.y)^2+(y.x-y.y)^2+(z.x-z.y)^2))
        df_test<-df_test%>%filter(length<12)
        df_test<-df_test%>%filter(alt.x!=alt.y)
        df_test<-df_test%>%select(eleno.x,  elety.x,  alt.x,    resid.x,  resno.x,x.x,y.x,z.x,
                                  eleno.y,  elety.y,  alt.y,    resid.y,  resno.y,x.y,y.y,z.y, length)
        df_test<-df_test%>%mutate(elec=NA)
        df_test<-df_test%>%mutate(vdw=NA)
        df_test<-df_test%>%mutate(nonb=NA)
        df_test<-df_test%>%mutate(Total=NA)

        
        df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
        df_tcl[1,1]<-paste0('cd ', part_start,"fin_structure/",v_ligand[i],"/\n\npackage require namdenergy")
        df_tcl[1,2]<-paste0('mol new {full_structure/',df_complex_TEMP$complex_name[j],'.psf} type {psf}')
        df_tcl[1,3]<-paste0('mol addfile {full_structure/',df_complex_TEMP$complex_name[j],'.pdb} type {pdb}')
        for (q in 1:nrow(df_test)) {
          df_tcl[(q+1),1]<-paste0('set sel2 [atomselect top "resid ',df_test$resno.x[q], ' and index ',df_test$eleno.x[q], '"]')
          df_tcl[(q+1),2]<-paste0('set sel1 [atomselect top "resid ',df_test$resno.y[q], ' and index ',df_test$eleno.y[q], '"]')
          df_tcl[(q+1),3]<-paste0('\nnamdenergy -sel $sel1 $sel2 -vdw -elec -nonb -cutoff 12 -skip 0 -ofile ',df_complex_TEMP$complex_name[j],'/bond_',
                                  df_test$resno.x[q],"-",df_test$eleno.x[q],"-",df_test$elety.x[q],"_",
                                  df_test$resno.y[q],"-",df_test$eleno.x[q],"-",df_test$elety.y[q],
                                  '.txt -switch 10 -exe /home/nastia/NAMD_2.14_Linux-x86_64-multicore/namd2 -par ',part_start,
                                  'start/toppar/par_all36_carb.prm -par ',part_start,'start/toppar/par_all36_cgenff.prm -par ',part_start,
                                  'start/toppar/par_all36_lipid.prm -par ',part_start,'start/toppar/par_all36m_prot.prm -par ',part_start,
                                  'start/toppar/par_all36_na.prm -par ',part_start,'start/toppar/par_all36_prot.prm -par ',part_start,
                                  'start/toppar/toppar_water_ions_namd.str')
          
         }
        df_tcl[(q+2),1]<-paste0('mol delete all \n \n \n exit now')
        write.table(df_tcl,paste0('fin_structure/tcl/energy_full_structure_',df_complex_TEMP$complex_name[j],'.tcl'),na = "\n",col.names = F,row.names = F,quote = F,sep = "\n")
        system(command = paste0("vmd -dispdev text -e ",part_start,'fin_structure/tcl/energy_full_structure_',df_complex_TEMP$complex_name[j],'.tcl'),ignore.stdout=T,wait = T) 
        q<-1
        for (q in 1:nrow(df_test)) {
          if(file.exists(paste0("fin_structure/",v_ligand[i],"/",df_complex_TEMP$complex_name[j],'/bond_',df_test$resno.x[q],"-",df_test$eleno.x[q],"-",df_test$elety.x[q],"_",
                                df_test$resno.y[q],"-",df_test$eleno.x[q],"-",df_test$elety.y[q],'.txt'))){
            df_energy<-read.table(paste0("fin_structure/",v_ligand[i],"/",df_complex_TEMP$complex_name[j],'/bond_',df_test$resno.x[q],"-",df_test$eleno.x[q],"-",df_test$elety.x[q],"_",
                                         df_test$resno.y[q],"-",df_test$eleno.x[q],"-",df_test$elety.y[q],'.txt'), sep="", header=T, na.strings ="", stringsAsFactors= F)
            df_test$elec[q]<-df_energy$Elec[1]
            df_test$vdw[q]<-df_energy$VdW[1]
            df_test$nonb[q]<-df_energy$Nonbond[1]
            df_test$Total[q]<-df_energy$Total[1]
          }
        }
        df_test_E<-df_test%>%filter(elec<(-0.1))
        df_test_E<-df_test_E%>%mutate(bond_type="Elec")
        df_test_W<-df_test%>%filter(vdw>(0.1))
        df_test_W<-df_test_W%>%mutate(bond_type="VdW")
        df_testa<-rbind(df_test_E,df_test_W)
        df_testa<-df_testa%>%mutate(sort_string=paste(eleno.y,resno.x))
        df_testa<-df_testa%>%group_by(sort_string)%>%mutate(min_Total=min(Total))
        df_testa<-df_testa%>%filter(Total==min_Total)
        df_testa<-ungroup(df_testa)
        df_testa<-df_testa%>%mutate(sort_string=paste(eleno.x,resno.x))
        df_testa<-df_testa%>%group_by(sort_string)%>%mutate(min_Total=min(Total))
        df_testa<-df_testa%>%filter(Total==min_Total)
        df_testa<-ungroup(df_testa)
        df_testa<-df_testa%>%select(eleno.x,elety.x, resid.x, resno.x, x.x,  y.x,  z.x,  
                                    eleno.y,elety.y, resid.y, resno.y, x.y,  y.y,  z.y,
                                    length, elec,    vdw,     nonb, Total, bond_type)
        write.csv(df_testa,paste0('fin_structure/interactions_atom/',df_complex_TEMP$complex_name[j],".csv"),row.names = F)
      }
    }
  }
}