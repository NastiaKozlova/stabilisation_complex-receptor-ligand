part_start = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
library(ggplot2)
library(cowplot)
name<-strsplit(part_start,split = "/",fixed = T)[[1]]
based_name<-name[length(name)]
num_model<-3
name<-based_name
setwd(part_start)
part_start<-paste0(part_start,'MD/stabilisation/')
setwd(part_start)
j<-1
for (j in 1:length(name)) {
  part_name<-paste0(part_start,"din/")
  setwd(part_name)
  if(file.exists(paste0("RMSF/",name[j],".txt"))){
    df_RMSF <- read.table(paste0("RMSF/",name[j],".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
    colnames(df_RMSF)<-"RMSF"
    df_RMSF <- df_RMSF%>% mutate(Resid=c(1:nrow(df_RMSF)))
    
    df_RMSD <- read.table(paste0("RMSD/",name[j],".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
    colnames(df_RMSD)<-c("frame","RMSD")
    #  df_RMSD<-df_RMSD%>% mutate(Time=Time/100) 
    df_second<-read.csv(paste0("pdb_second/Second_structure_",name[j],".csv"), stringsAsFactors= F)
    df_SASA <- read.table(paste0("SASA/",name[j],".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
    colnames(df_SASA)<-c("frame","protein")
    df_energy<-read.table(paste0("Energy/protein_",name[j],".txt"), sep="", header=T, na.strings ="", stringsAsFactors= F)
    df_SASA<-df_SASA%>%mutate(frame = 0:(nrow(df_SASA)-1))
    df_energy<-df_energy%>%mutate(frame = 0:(nrow(df_energy)-1))
    df_ramachadran<-read.csv(paste0("pdb_second/",name[j],"_time_Ramachadran.csv"),stringsAsFactors = F)
    df_data<-left_join(df_ramachadran,df_energy,by=c("number"="frame"))
    df_data<-left_join(df_data,df_SASA,by=c("number"="frame"))
    df_data<-left_join(df_data,df_RMSD,by=c("number"="frame"))
    df_data<-df_data%>%mutate(time=number/100)
    df_second<-df_second%>%mutate(level_min=level_min/100)
    df_second<-df_second%>%mutate(level_max=level_max/100)
    v_RMSD<-df_data$RMSD[!is.na(df_data$RMSD)]
    p_rmsd<-ggplot(data = df_data)+
      ggtitle(paste0("RMSD"))+
      labs(y = "RMSD (A)", x = "Time (ns)")+
      labs(y = "RMSD (A)", x = "Время, нс")+
      geom_line(aes(x = time, y = RMSD))+
      geom_hline(yintercept = median(v_RMSD))+
      #    scale_x_continuous(breaks = test_10, labels =  test_10)+
      scale_y_continuous(limits = c(min(v_RMSD),max(v_RMSD)))+
      theme_bw()+coord_flip()
    v_protein<-df_data$protein[!is.na(df_data$protein)]
    p_sasa<-ggplot(data = df_data)+
      labs(title=paste("Solvent accessible surface"), x = "Time(ns)", y = "Solvent accessible surface (A^2)") +
      labs(title=paste("Площадь доступная для растрорителя"),x =  "Время, нс", y = "Площадь доступная для растрорителя (A^2)") +
      geom_line(aes(x = time,y=protein))+
      geom_hline(yintercept = median(v_protein))+
      theme_bw()+coord_flip()
    v_Total<-df_data$Total[!is.na(df_data$Total)]
    p_Total<-ggplot(data = df_data)+
      labs(title=paste("Total energy"), x = "Time(ns)", y = "Total energy (kcal/mol)") +
      labs(title=paste("Полная энергия рецептора"), x = "Время, нс", y = "Полная энергия рецептораб ккал/моль") +
      geom_line(aes(x = time,y=Total))+
      geom_hline(yintercept = median(v_Total))+
      #    scale_x_continuous(breaks = test_10, labels =  test_10)+
      #    scale_y_continuous(limits = c(min(df_energy$Total),max(df_energy$Total)))+
      theme_bw()+coord_flip()
    v_ramachadran<-df_data$ramachadran[!is.na(df_data$ramachadran)]
    p_ramachadran<-ggplot(data = df_data)+
#      labs(title=paste("Ramachadran"),
#           x = "Time(ns)", y = "Ramachadran") +
      labs(title=paste("Кол-во аминокислот в запрещенных зонах"),x = "Время, нс", y = "Кол-во аминокислот в запрещенных зонах") +
      geom_line(aes(x = time,y=ramachadran))+
      geom_hline(yintercept = median(v_ramachadran))+
      #    scale_x_continuous(breaks = test_10, labels =  test_10)+
      #    scale_y_continuous(limits = c(min(df_ramachadran$ramachadran),max(df_ramachadran$ramachadran)))+
      theme_bw()+coord_flip()
    
    p_ramachadran_histo<-ggplot(data = df_data)+
      labs(#title=paste("Ramachadran"),
           x = "Time(ns)", y = "Ramachadran") +
      geom_freqpoly(aes(x = ramachadran),bins=(max(v_ramachadran)-min(v_ramachadran)))+
      geom_vline(xintercept = median(v_ramachadran))+
      geom_text(aes(x = median(v_ramachadran),y=100,label= round(median(v_ramachadran),digits = 0)))+
      #     scale_x_continuous(breaks = test_10, labels =  test_10)+
      scale_x_continuous(limits = c(min(v_ramachadran),max(v_ramachadran)))+
      theme_bw()
    
    p_rmsd_histo<-ggplot(data = df_data)+
#      ggtitle(paste0("RMSD"))+
      labs(x = "RMSD (A)")+
      #      geom_freqpoly(aes(x = RMSD),bins=(max(df_RMSD$RMSD)-min(df_RMSD$RMSD)))+
      geom_freqpoly(aes(x = RMSD))+
      geom_vline(xintercept = median(v_RMSD))+
      geom_text(aes(x = median(v_RMSD),y=100,label= round(median(v_RMSD),digits = 0)))+
      scale_x_continuous(limits = c(min(v_RMSD),max(v_RMSD)))+
      theme_bw()
    p_sasa_histo<-ggplot(data = df_data)+
      labs(#title=paste("Solvent accessible surface"),
           x =  "Solvent accessible surface (A^2)") +
      #      geom_freqpoly(aes(x = protein),bins=(max(df_SASA$protein)-min(df_SASA$protein)))+
      geom_freqpoly(aes(x = protein))+
      geom_vline(xintercept = median(v_protein))+
      geom_text(aes(x = median(v_protein),y=100,label= round(median(v_protein),digits = 0)))+
      scale_x_continuous(limits = c(min(v_protein),max(v_protein)))+
      theme_bw()
    p_Total_histo<-ggplot(data = df_data)+
      labs(#title=paste("Total energy"),
#           x =  "Total energy (kcal/mol)" +
      x =  "Полная энергия, ккал/моль") +
      geom_freqpoly(aes(x = Total))+
      geom_vline(xintercept = median(v_Total))+
      geom_text(aes(x = median(v_Total),y=100,label= round(median(v_Total),digits = 0)))+
      scale_x_continuous(limits = c(min(v_Total),max(v_Total)))+
      theme_bw()
    
    
    
    p_second<-ggplot(data = df_second)+
      ggtitle(paste0("Second structure"))+
      labs(x = "Number of aminoasids", y = "Time, ns")+       labs(x = "Номер аминокислоты", y = "Время, нс")+
      geom_rect(aes(xmin = start, ymin = level_min, xmax= finish, ymax = level_max, colour = type,fill=type))+
      scale_color_grey()+ scale_fill_grey()+
      theme_bw()+ guides(color = "none") + guides(fill = "none")
    
    p_rmsf<-ggplot(data = df_RMSF)+
      ggtitle(paste0("RMSF"))+
      labs(x = "Number of aminoasids", y = "RMSF (A)")+labs(x = "Номер аминокислоты", y = "RMSF (A)")+
      geom_line(aes(x = Resid, y = RMSF))+
      theme_bw()
    p_all<-plot_grid(p_sasa,       p_ramachadran,       p_Total,       p_rmsd,       p_second,
                     p_sasa_histo, p_ramachadran_histo, p_Total_histo, p_rmsd_histo, p_rmsf, nrow=2,rel_heights = c(4,1),rel_widths = c(1,1,1,1,4.5),align = "hv")
    ggsave(p_all,filename = paste0("RMSD_RMSF_Second_str_",name[j],".png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 
    write.csv(df_data,"df_select_number.csv",row.names = F)
    m_SASA<-quantile(df_data$protein,probs = 0.025)
    m_Total<-quantile(df_data$Total,probs = 0.025)
    df_data<-df_data%>%filter(ramachadran==min(df_data$ramachadran))
    df_data<-df_data%>%filter(Total<m_Total)
    df_data<-df_data%>%filter(protein<m_SASA)
    df_RMSD<-df_data%>%select(number,frame_number)
    df_RMSD<-df_RMSD%>%mutate(RMSD=NA)
    df_RMSD_all<-full_join(df_RMSD,df_RMSD,by="RMSD")
    df_RMSD_all<-df_RMSD_all%>%filter(number.x<number.y)
    for (i in 1:nrow(df_RMSD_all)) {
      pdb_1<-read.pdb(paste0("pdb_second/",name[j],"/",df_RMSD_all$frame_number.x[i],".pdb"))
      pdb_2<-read.pdb(paste0("pdb_second/",name[j],"/",df_RMSD_all$frame_number.y[i],".pdb"))
      df_RMSD_all$RMSD[i]<-rmsd(pdb_1,pdb_2,fit = T)
    }
    write.csv(df_RMSD_all,"df_RMSD_all.csv",row.names =  F)
    
    df_RMSD_sep<-read.csv("df_RMSD_all.csv",stringsAsFactors = F)
    df_RMSD_control<-df_RMSD
    if (!file.exists(paste0('grop'))) {dir.create(paste0("grop"))}
    if (!file.exists(paste0('fin_structure'))) {dir.create(paste0("fin_structure"))}
    df_RMSD_sep<-df_RMSD_sep%>%filter(RMSD<5)
    for (i in 1:nrow(df_RMSD)) {
      if (!is.na(df_RMSD$number[i])) {
        df_RMSD_sep_test<-df_RMSD_sep%>%filter(number.x==df_RMSD$number[i])
        if (nrow(df_RMSD_sep_test)>(nrow(df_RMSD_control)/10)) {
          if (!file.exists(paste0('fin_structure/',i))) {dir.create(paste0("fin_structure/",i))}
          df_RMSD_sep_test_add<-data.frame(matrix(ncol=ncol(df_RMSD_sep_test),nrow = 1))
          colnames(df_RMSD_sep_test_add)<-colnames(df_RMSD_sep_test)
          df_RMSD_sep_test_add$number.x<-df_RMSD$number[i]
          df_RMSD_sep_test_add$number.y<-df_RMSD$number[i]
          df_RMSD_sep_test_add$frame_number.x<-df_RMSD$frame_number[i]
          df_RMSD_sep_test_add$frame_number.y<-df_RMSD$frame_number[i]
          df_RMSD_sep_test_add$number.y<-df_RMSD$number[i]
          df_RMSD_sep_test_add$RMSD<-0
          df_RMSD_sep_test<-rbind(df_RMSD_sep_test,df_RMSD_sep_test_add)
          pdb_name<-unique(df_RMSD_sep_test$frame_number.x)
          pdb_name<-pdb_name[!is.na(pdb_name)]
          write.csv(df_RMSD_sep_test,paste0("grop/",i,".csv"),row.names = F) 
          pdb<-read.pdb(paste0("pdb_second/",name[j],"/",pdb_name,".pdb"))
          write.pdb(pdb,paste0(part_start,"start/ligand_start/",name[j],"_",i,".pdb"))
        }
      }
      df_RMSD_sep$number.x[df_RMSD_sep$number.x%in%df_RMSD_sep_test$number.y]<-NA
      df_RMSD_sep$number.x[df_RMSD_sep$number.y%in%df_RMSD_sep_test$number.y]<-NA
      df_RMSD_sep<-df_RMSD_sep%>%filter(!is.na(number.x))
      df_RMSD$number[df_RMSD$number%in%df_RMSD_sep_test$number.y]<-NA
    }
  }
}