part_start <- commandArgs(trailingOnly=TRUE)
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
      
      if(!dir.exists(paste0("fin_structure/",v_ligand[i]))){dir.create((paste0("fin_structure/",v_ligand[i])))}
      for (j in 1:nrow(df_complex_TEMP)) {
        df_interactions<-read.csv(paste0(df_complex_TEMP$complex_name[j],"/MD/stabilisation/din/intercations.csv"),stringsAsFactors = F)
        df_interactions<-df_interactions%>%filter(number<1000)
        df_interactions<-df_interactions%>%filter(Total==min(df_interactions$Total))
        df_tcl[(j+1),1]<-paste0('mol new {prepared_structures/complex_structure/',df_complex_TEMP$complex_name[1],'.psf} type {psf}')
        df_tcl[(j+1),2]<-paste0('mol addfile {',df_complex_TEMP$complex_name[j],'/MD/stabilisation/din/pdb_second_complex/',df_complex_TEMP$complex_name[j],'/frame_',df_interactions$number[1],'.pdb} type {pdb}') 
        df_tcl[(j+1),3]<-paste0('[atomselect top all] writepdb fin_structure/',v_ligand[i],"/",df_complex_TEMP$complex_name[j],'.pdb\n')
        df_tcl[(j+1),4]<-paste0('[atomselect top all] writepsf fin_structure/',v_ligand[i],"/",df_complex_TEMP$complex_name[j],'.psf\n')
        df_tcl[(j+1),5]<-paste0('mol delete all \n \n \n ')
      }
      df_tcl[(j+2),1]<-paste0('mol delete all \n \n \n exit now')
      write.table(df_tcl,paste0('fin_structure/tcl/',v_ligand[i],'.tcl'),na = "\n",col.names = F,row.names = F,quote = F,sep = "\n")
      system(command = paste0("vmd -dispdev text -e ",part_start,'fin_structure/tcl/',v_ligand[i],'.tcl'),ignore.stdout=T,wait = T) 
    
      for (j in 1:nrow(df_complex_TEMP)) {
        df_interactions<-read.csv(paste0(df_complex_TEMP$complex_name[j],"/MD/stabilisation/din/intercations.csv"),stringsAsFactors = F)
        df_interactions<-df_interactions%>%filter(number<1000)
        df_interactions<-df_interactions%>%filter(Total==min(df_interactions$Total))
        df_interaction<-read.csv(paste0(df_complex_TEMP$complex_name[j],'/MD/stabilisation/din/interaction/',df_complex_TEMP$complex_name[j],'/frame_',df_interactions$number[1],'.csv'),stringsAsFactors = F)
        print(paste0(df_complex_TEMP$center[j]," resid ",paste0(df_interaction$resid,collapse = " ")))
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
        df_test<-df_test%>%filter(length<5)
        df_test<-df_test%>%filter(alt.x!=alt.y)
        df_test<-df_test%>%select(eleno.x,  elety.x,  alt.x,    resid.x,  resno.x,x.x,y.x,z.x,
                                  eleno.y,  elety.y,  alt.y,    resid.y,  resno.y,x.y,y.y,z.y, length)
        write.csv(df_test,paste0('fin_structure/interactions/',v_ligand[i],"_",df_complex_TEMP$complex_name[j],".csv"),row.names = F)
      }
      df_interactions_start<-read.csv(paste0(df_complex_TEMP$complex_name[1],"/MD/stabilisation/din/intercations.csv"),stringsAsFactors = F)
      df_interactions_start<-df_interactions_start%>%mutate(complex_name=df_complex_TEMP$complex_name[1])
      if(nrow(df_complex_TEMP)>1){
        for (j in 2:nrow(df_complex_TEMP)) {
          df_interactions<-read.csv(paste0(df_complex_TEMP$complex_name[j],"/MD/stabilisation/din/intercations.csv"),stringsAsFactors = F)
          df_interactions<-df_interactions%>%mutate(complex_name=df_complex_TEMP$complex_name[j])
          df_interactions_start<-rbind(df_interactions_start,df_interactions)
        }
      }
      df_interactions_start<-df_interactions_start%>%group_by(complex_name)%>%mutate(total_median=median(Total))
      
      df_interactions_start<-left_join(df_interactions_start,df_complex_TEMP,by = "complex_name")
      df_interactions_start<-df_interactions_start%>%filter(Total!=0)
      df_interactions_start<-df_interactions_start%>%mutate(center_new=NA)
      df_interactions_start<-df_interactions_start%>%mutate(ligand_new=NA)
      df_interactions_start$center_new[df_interactions_start$center=="perefer_anion_cite"]<-"1. Периферийный анионный сайт"
      df_interactions_start$center_new[df_interactions_start$center=="acyl_poket"]<-"4. Ацильный карман"
      df_interactions_start$center_new[df_interactions_start$center=="anion_cite"]<-"3. Анионный сайт"
      df_interactions_start$center_new[df_interactions_start$center=="hydroxyanion_hole"]<-"2. Оксианионная дыра"
      df_interactions_start$center_new[df_interactions_start$center=="cat_trio"]<-"5. Каталитическая триада"

      df_interactions_start$ligand_new[df_interactions_start$ligand=="ACh"]<-"АХ"
      df_interactions_start$ligand_new[df_interactions_start$ligand=="ATC"]<-"АТХ"
      df_interactions_start$ligand_new[df_interactions_start$ligand=="atma"]<-"АТМА"
      df_interactions_start$ligand_new[df_interactions_start$ligand=="BCh"]<-"БХ"
      df_interactions_start$ligand_new[df_interactions_start$ligand=="BTC"]<-"БТХ"
      df_interactions_start$ligand_new[df_interactions_start$ligand=="BzCh"]<-"БзХ"
      df_interactions_start$ligand_new[df_interactions_start$ligand=="BzTC"]<-"БзТХ"
      df_interactions_start$ligand_new[df_interactions_start$ligand=="choline"]<-"холин"
      
      p<-ggplot(data = df_interactions_start)+
        labs(x="Аффинность, ккал/моль",y="")+
        geom_freqpoly(aes(x=Total),binwidth=10)+
        geom_text(aes(x=total_median,y=10,label=total_median))+
        facet_grid(.~center_new)+theme_bw()
      ggsave(p,filename = paste0("fin_structure/plot/fin_structure_energy_",v_ligand[i],".png"), width = 60, height = 15, units = c("cm"), dpi = 200 ) 
    }
  }
}