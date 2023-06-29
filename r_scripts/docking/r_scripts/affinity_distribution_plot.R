part_analysis <- commandArgs(trailingOnly=TRUE)
#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(ggpmisc)
#v_rmsd<-10
#v_group_size<-24

#setwd(part_analysis)
#df_all<-read.csv(paste0(part_analysis,"df_all.csv"),stringsAsFactors = F)
#df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))

part<-paste0(part_analysis,"din/")
setwd(part)
df_structure_RMSD<-read.csv("df_merge_structure_log.csv",stringsAsFactors = F)
df_structure_RMSD<-df_structure_RMSD[df_structure_RMSD$ligand%in%c("DMPE","DPPE",
                                                                   "DPPG","DYPE",
                                                                   "DYPG","PYPE"),]
# Ненормальное распределение
ggboxplot(df_structure_RMSD, x = "receptor_ligand", y = "affinity")

ggqqplot(df_structure_RMSD, "affinity", facet.by = "receptor_ligand")

kruskal <- df_structure_RMSD %>% kruskal_test(affinity ~ receptor_ligand)
# Не достоверно!

pwc2 <- df_structure_RMSD %>% 
   dunn_test(affinity ~ receptor_ligand, p.adjust.method = "bonferroni") 

pwc2 <- pwc2 %>% add_xy_position(x = "receptor_ligand")
pwc2<-pwc2%>%filter(p.adj.signif!="ns")
pwc2<-pwc2%>%filter(p.adj.signif!="*")
p<-ggboxplot(df_structure_RMSD, x = "receptor_ligand", y = "affinity") +
   stat_pvalue_manual(pwc2) +
   theme_bw()+
   labs(
      subtitle = get_test_label(kruskal, detailed = TRUE),
      caption = get_pwc_label(pwc2)
   )

ggsave(p,file="boxplot_ligand_affility.png", width = 24, height = 15, units = c("cm"), dpi = 200 )
colnames(pwc2)
pwc2<-pwc2%>%select(group1,group2,p.adj,p.adj.signif)
p<-ggplot(data=df_structure_RMSD, aes(colour = receptor_ligand, x = affinity)) +
   geom_density()+theme_bw()+
   theme(legend.position = "bottom")
p_table<-p+annotate(geom = "table",
                  x = 0,
                  y = 0.5,
                  label = list(pwc2))
#pwc2
ggsave(p_table,file="density_ligand_affility.png", width = 24, height = 15, units = c("cm"), dpi = 200 )
