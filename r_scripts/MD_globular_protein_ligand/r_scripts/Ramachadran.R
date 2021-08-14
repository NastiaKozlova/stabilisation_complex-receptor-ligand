part_start = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(bio3d)
library(ggplot2)
setwd(part_start)
name<-strsplit(part_start,split = "/",fixed = T)[[1]]
part_rscrpts<-paste0(name[1:(length(name)-1)],collapse = "/")
based_name<-name[length(name)]
num_model<-3
name<-1
part_start<-paste0(part_start,'MD/stabilisation/')

make_df_ramachadran<-function(EMD2){
  beg<-EMD2$atom$resno[1]
  fin<-EMD2$atom$resno[length(EMD2$atom$resno)]
  tor<-torsion.pdb(EMD2)
  df_tor<-data.frame(matrix(nrow = length(tor$phi),ncol=2))
  colnames(df_tor)<-c("phi","psi")
  df_tor[1:length(tor$phi),1]<-tor$phi
  df_tor[1:length(tor$psi),2]<-tor$psi
  df_tor$seq<-NA
  df_tor$seq[1:length(EMD2$atom$resid[atom.select(EMD2,"calpha")$atom])]<-EMD2$atom$resid[atom.select(EMD2,"calpha")$atom]
  df_tor$number<-NA
  df_tor$number<-c(beg:fin)
  df_tor<- df_tor%>% filter(seq!="GLY")
  df_tor<-df_tor%>%mutate(amino=paste(number,seq))
  return(df_tor)
}

filter_df_ramachadran<-function(df_rama,df_rama_start){
  df_rama<-df_rama%>%filter(!is.na(phi))
  df_rama<-df_rama%>%filter(!is.na(psi))
  
  df_rama<-df_rama%>%mutate(favorite=0)
  o<-0
  for (o in 0:(ncol(df_rama_start)/2-1)) {
    
    df_rama_test<-df_rama_start[,(2*o+1):(2*o+2)]
    colnames(df_rama_test)<-c("phi","psi")
    df_rama_test<-df_rama_test%>%filter(!is.na(phi))
    df_rama_test<-df_rama_test%>%filter(!is.na(psi))
    
    for (a in 1:nrow(df_rama)) {
      x0<-df_rama$phi[a]
      y0<-df_rama$psi[a]
      p<-0
      
      for (j in 2:nrow(df_rama_test)) {
        q<- ((df_rama_test$phi[j-1]-x0)*(df_rama_test$phi[j]-x0)+(df_rama_test$psi[j-1]-y0)*(df_rama_test$psi[j]-y0))
        w<- (sqrt((df_rama_test$phi[j-1]-x0)^2+(df_rama_test$psi[j-1]-y0)^2)*sqrt((df_rama_test$phi[j]-x0)^2+(df_rama_test$psi[j]-y0)^2))
        g<-q/w
        if (abs(g)>1) {g<-(1) }
        e<- round(acos(g),digits = 3)
        r<-(-(df_rama_test$phi[j-1]-x0)*(df_rama_test$psi[j]-y0)+(df_rama_test$psi[j-1]-y0)*(df_rama_test$phi[j]-x0))
        
        if (r<0) {p<-p-e  } else{ p<-p+e}
        
        if (abs(p)>df_rama$favorite[a]) {
          df_rama$favorite[a]<-abs(p)
        }
      }
    } 
  }
  df_rama_filter<-df_rama%>%filter(favorite<6)
  return(df_rama_filter)
}

part<-part_start
p<-1
for (p in 1:length(based_name)) {
  part_start<-paste0(part,"/din/pdb_second/")
  setwd(part_start)
  if (file.exists(paste0(based_name[p],"/frame_",0,".pdb"))){ 
    frame_number<-length(list.files(path = paste0(based_name[p])))
    df_topology<-data.frame(matrix(nrow = frame_number, ncol=3))
    colnames(df_topology)<-c("number","frame_number","ramachadran")
    df_topology$number<-0:(frame_number-1)
    df_topology<-df_topology%>%mutate(frame_number=paste0("frame_",number))
    if(!dir.exists(paste0(based_name[p],"_rama/"))){(dir.create(paste0(based_name[p],"_rama/")))}
    for (i in 1:nrow(df_topology)) {
      protein<-read.pdb(paste0(based_name[p],"/",df_topology$frame_number[i],".pdb"))
      protein_atom<-atom.select(protein,"protein")
      EMD2<-trim.pdb(protein,protein_atom)
      df_rama<-make_df_ramachadran(EMD2)
      if(!dir.exists(paste0(based_name[p],"_rama/")))    { dir.create(paste0(based_name[p],"_rama/"))}
      write.csv(df_rama,file = paste0(based_name[p],"_rama/",df_topology$frame_number[i],".csv"),row.names = F)
    }
    df_rama_start<-read.csv(paste0(part_rscrpts,"/r_scripts/rama4.csv"),stringsAsFactors = F)
    for (i in 1:nrow(df_topology)) {
      df_rama<-read.csv(file = paste0(based_name[p],"_rama/",df_topology$frame_number[i],".csv"),stringsAsFactors =  F)
      if(!dir.exists(paste0(based_name[p],"_ramachad/")))    { dir.create(paste0(based_name[p],"_ramachad/"))}
      df_rama_filter<-filter_df_ramachadran(df_rama = df_rama,df_rama_start = df_rama_start)
      write.csv(df_rama_filter,paste0(based_name[p],"_ramachad/",df_topology$frame_number[i],".csv"),row.names = F)
    }
    for (i in 1:nrow(df_topology)) {
      df_rama_filter<-read.csv(paste0(based_name[p],"_ramachad/",df_topology$frame_number[i],".csv"),stringsAsFactors = F)
      df_topology$ramachadran[i] <- nrow(df_rama_filter)
    }
    df_topology<-df_topology%>%mutate(time=number/1000)
    p1<-ggplot(data = df_topology)+
      labs(y="Nimber of aminoasid", x="Time, ns")+
      geom_line(aes(x=time,y=ramachadran))+theme_bw()
    ggsave(plot = p1, filename = paste0(based_name[p],"_time_Ramachadran.png"),  width = 20, height = 20, units = c("cm"), dpi = 300 )
    write.csv(df_topology,paste0(based_name[p],"_time_Ramachadran.csv"),row.names = F)
    df_topology<-df_topology%>%group_by(ramachadran)%>%mutate(x_num=n())
    df_topology<-ungroup(df_topology)
    df_topology<-df_topology%>%select(ramachadran, x_num)
    df_topology<-unique(df_topology)
    df_topology<-df_topology%>%mutate(x_num=x_num/frame_number*100)
    write.csv(df_topology,paste0(based_name[p],"_Ramachadran.csv"),row.names = F)
  }
  if (file.exists(paste0(based_name[p],"_Ramachadran.csv"))){
    df_topology<-read.csv(paste0(based_name[p],"_Ramachadran.csv"),stringsAsFactors =  F)
    df_topology<-df_topology%>%arrange(ramachadran)
    df_topology<-df_topology%>%mutate(persent=x_num)
    for (i in 2:nrow(df_topology)) {
      df_topology$persent[i]<-df_topology$persent[i-1]+df_topology$x_num[i]
    }
    df_topology<-df_topology%>%mutate(persent=round(persent,digits = 2))
    df_topology$persent[df_topology$persent>1]<-round(df_topology$persent[df_topology$persent>1],digits = 1)
    p1<-ggplot(data = df_topology)+
      labs(x="Quantity of unfavorite aminoasids",y="persent of time, %")+
      geom_point(aes(x=ramachadran,y=x_num))+
      geom_line(aes(x=ramachadran,y=x_num))+theme_bw()+
      geom_text(aes(x=ramachadran,y=x_num+0.25,label=persent))+
      scale_y_continuous(breaks = c(min(round(df_topology$x_num,digits = 0)):max(round(df_topology$x_num,digits = 0))),
                         labels = c(min(round(df_topology$x_num,digits = 0)):max(round(df_topology$x_num,digits = 0))))+
      scale_x_continuous(breaks = c(min(df_topology$ramachadran):max(df_topology$ramachadran)),
                         labels = c(min(df_topology$ramachadran):max(df_topology$ramachadran)))
    ggsave(plot = p1, filename = paste0(based_name[p],"_Ramachadran.png"),  width = 20, height = 15, units = c("cm"), dpi = 300 )
    df_topology<-read.csv(paste0(based_name[p],"_time_Ramachadran.csv"),stringsAsFactors = F)
    df_rama<-read.csv(paste0(based_name[p],"_ramachad/",df_topology$frame_number[1],".csv"),stringsAsFactors =  F)
    for (i in 2:nrow(df_topology)) {
      df_rama_filter<-read.csv(paste0(based_name[p],"_ramachad/",df_topology$frame_number[i],".csv"),stringsAsFactors = F)
      df_rama<-rbind(df_rama,df_rama_filter)
    }
    df_rama<-df_rama%>%group_by(amino)%>%mutate(x_num=n())
    df_rama<-ungroup(df_rama)
    df_rama<-df_rama%>%select(seq, number ,amino, x_num)
    df_rama<-unique(df_rama) 
    df_rama<-df_rama%>%mutate(x_num=x_num/frame_number*100)
    df_rama_filter<-df_rama
    p1<-ggplot(data = df_rama)+
      labs(x="Nimber of aminoasid", y="persent of time, %")+
      geom_point(aes(x=number,y=x_num))+theme_bw()+
      geom_text(aes(x=number,y=x_num,label=amino),data = df_rama_filter)
    write.csv(df_rama,paste0(based_name[p],"_ramachad.csv"),row.names = F)
    ggsave(plot = p1, filename = paste0(based_name[p],"_Ramachadran_histogram.png"),  width = 40, height = 20, units = c("cm"), dpi = 300 )
  }
}
