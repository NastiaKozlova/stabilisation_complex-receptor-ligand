part_name = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(bio3d)
library(ggplot2)
library(cowplot)
#part_name<-part_start
setwd(part_name)
name<-list.files("start/structure/")

a<-c()
for(i in 1:length(name)){
  b<-strsplit(name[i],split = ".",fixed = T)[[1]][1]
  a<-c(a,b)
}
name<-a
j<-1
i<-1
for (j in 1:length(name)) {
  part<-paste0(part_name,name[j],"/MD/stabilisation/din/")
  if(dir.exists(part)){
    setwd(part)
    if(!dir.exists("fin_structure")){dir.create("fin_structure")}
    df_RMSF <- read.table(paste0("RMSF/",name[j],".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
    colnames(df_RMSF)<-"RMSF"
    df_RMSF <- df_RMSF%>% mutate(Resid=c(1:nrow(df_RMSF)))
    
    df_RMSD <- read.table(paste0("RMSD/",name[j],".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
    colnames(df_RMSD)<-c("frame","RMSD")
    #  df_RMSD<-df_RMSD%>% mutate(Time=Time/10) 
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
    df_data<-df_data%>%mutate(time=number/10)
    df_second<-df_second%>%mutate(level_min=level_min/10)
    df_second<-df_second%>%mutate(level_max=level_max/10)
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
      labs(title=paste("Solvent accessible surface"),
           x = "Time(ns)", y = "Solvent accessible surface (A^2)") +
      labs(title=paste("Площадь доступ. для раств-ля"),
           x="Время, нс",y =  "Площадь доступ. для раств-ля,A^2") +
      geom_line(aes(x = time,y=protein))+
      geom_hline(yintercept = median(v_protein))+
      #    scale_x_continuous(breaks = test_10, labels =  test_10)+
      #   scale_y_continuous(limits = c(min(df_SASA$protein),max(df_SASA$protein)))+
      theme_bw()+coord_flip()
    v_Total<-df_data$Total[!is.na(df_data$Total)]
    p_Total<-ggplot(data = df_data)+
      labs(title=paste("Total energy"),
           x = "Time(ns)", y = "Total energy (kcal/mol)") +
      labs(title=paste("Полная энергия"),
           x ="Время, нс",y = "Полная энергия, ккал/моль") +
      geom_line(aes(x = time,y=Total))+
      geom_hline(yintercept = median(v_Total))+
      #    scale_x_continuous(breaks = test_10, labels =  test_10)+
      #    scale_y_continuous(limits = c(min(df_energy$Total),max(df_energy$Total)))+
      theme_bw()+coord_flip()
    v_ramachadran<-df_data$ramachadran[!is.na(df_data$ramachadran)]
    p_ramachadran<-ggplot(data = df_data)+
      labs(title=paste("Ramachadran"),
           x = "Time(ns)", y = "Ramachadran") +
      labs(title=paste("Кол-во аминокслот в запрещенных зонах"),
           x = "Время, нс", y = "Кол-во аминокслот в запрещенных зонах") +
      geom_line(aes(x = time,y=ramachadran))+
      geom_hline(yintercept = median(v_ramachadran))+
      #    scale_x_continuous(breaks = test_10, labels =  test_10)+
      #    scale_y_continuous(limits = c(min(df_ramachadran$ramachadran),max(df_ramachadran$ramachadran)))+
      theme_bw()+coord_flip()
    
    p_ramachadran_histo<-ggplot(data = df_data)+
      labs(x = "Ramachadran",y="") +
      labs(x = "Кол-во аминокслот в запрещенных зонах") +
      geom_freqpoly(aes(x = ramachadran),bins=(max(v_ramachadran)-min(v_ramachadran)))+
      geom_vline(xintercept = median(v_ramachadran))+
      geom_text(aes(x = median(v_ramachadran),y=100,label= round(median(v_ramachadran),digits = 0)))+
      #     scale_x_continuous(breaks = test_10, labels =  test_10)+
      scale_x_continuous(limits = c(min(v_ramachadran),max(v_ramachadran)))+
      theme_bw()
    
    p_rmsd_histo<-ggplot(data = df_data)+
      labs(x = "RMSD (A)",y="")+
      #      geom_freqpoly(aes(x = RMSD),bins=(max(df_RMSD$RMSD)-min(df_RMSD$RMSD)))+
      geom_freqpoly(aes(x = RMSD))+
      geom_vline(xintercept = median(v_RMSD))+
      geom_text(aes(x = median(v_RMSD),y=100,label= round(median(v_RMSD),digits = 0)))+
      scale_x_continuous(limits = c(min(v_RMSD),max(v_RMSD)))+
      theme_bw()
    p_sasa_histo<-ggplot(data = df_data)+
      labs(x =  "Solvent accessible surface (A^2)") +
      labs(x =  "Площадь доступ. для раств-ля,A^2",y="")+
      geom_freqpoly(aes(x = protein))+
      geom_vline(xintercept = median(v_protein))+
      geom_text(aes(x = median(v_protein),y=100,label= round(median(v_protein),digits = 0)))+
      scale_x_continuous(limits = c(min(v_protein),max(v_protein)))+
      theme_bw()
    p_Total_histo<-ggplot(data = df_data)+
      labs(x =  "Total energy (kcal/mol)",y="") +
      labs(x =  "Полная энергия, ккал/моль",y="") +
      geom_freqpoly(aes(x = Total))+
      geom_vline(xintercept = median(v_Total))+
      geom_text(aes(x = median(v_Total),y=100,label= round(median(v_Total),digits = 0)))+
      scale_x_continuous(limits = c(min(v_Total),max(v_Total)))+
      theme_bw()
    
    
    
    p_second<-ggplot(data = df_second)+
      ggtitle(paste0("Second structure"))+
      ggtitle(paste0("Вторичная структура"))+
      labs(x = "Number of aminoasids", y = "Time, ns")+
      labs(x = "Номер аминокислоты", y = "Время, нс")+
      geom_rect(aes(xmin = start, ymin = level_min, xmax= finish, ymax = level_max, colour = type,fill=type))+
      scale_color_grey()+ scale_fill_grey()+
      theme_bw()+ guides(color = "none") + guides(fill = "none")
    
    p_rmsf<-ggplot(data = df_RMSF)+
      ggtitle(paste0("RMSF"))+
      labs(x = "Number of aminoasids", y = "RMSF (A)")+
      labs(x = "Номер аминокислоты", y = "RMSF (A)")+
      geom_line(aes(x = Resid, y = RMSF))+
      theme_bw()
    p_all<-plot_grid(p_sasa,       p_ramachadran,       p_Total,       p_rmsd,       p_second,
                     p_sasa_histo, p_ramachadran_histo, p_Total_histo, p_rmsd_histo, p_rmsf, nrow=2,rel_heights = c(4,1),rel_widths = c(1,1,1,1,4.5),align = "hv")
    ggsave(p_all,filename = paste0("RMSD_RMSF_Second_str_",name[j],".png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 
    
    write.csv(df_data,"df_select_number.csv",row.names = F)
    df_data<-read.csv("df_select_number.csv",stringsAsFactors = F)
    df_data<-df_data%>%filter(number>1)
    m_Ramachadran<-min(df_data$ramachadran)
    df_data<-df_data%>%filter(ramachadran==m_Ramachadran)
    df_data<-df_data%>%filter(Total==min(Total))
    df_data<-df_data%>%arrange(desc(time))
    pdb<-read.pdb(paste0("pdb_second/",name[j],"/",df_data$frame_number[1],".pdb"))
    write.pdb(pdb,paste0("fin_structure/",name[j],".pdb"))
  }
}