library(bio3d)
library(dplyr)
part_start = commandArgs(trailingOnly=TRUE)
setwd(part_start)
  if (!dir.exists(paste0(part_start,paste0(main,"/structure_prediction/pdb_fin")))){dir.create(paste0(part_start,paste0(main,"/structure_prediction/pdb_fin")))}
  seq<-read.fasta(paste0(main,".fasta"))
  fasta<-read.fasta(paste0(main,"/structure_prediction/seqs.fasta"))
  tost<-seqbind(seq,fasta)
  aln<-seqaln(tost,seqgroup = T, exefile="muscle") 
  df_seq<-as.data.frame(t((aln$ali)),stringsAsFactors = F)
  v_names<-colnames(df_seq)
  df_structure<-data.frame(matrix(ncol=2,nrow = length(v_names)))
  colnames(df_structure)<-c("structure","aminoasid_number")
  df_structure$structure<-v_names
  for (i in 1:length(v_names)) {
    if(length(strsplit(colnames(df_seq)[i],split = ":",fixed = T)[[1]])==1){
      start<-colnames(df_seq)[i]
      number<-i
    }
  }

  df_seqv<-df_seq%>%filter(df_seq[,number]!="-")
  i<-1
  for (i in 1:ncol(df_seqv)) {
#  for (i in 1:(nrow(df_seqv)-1)) {
    v_leng<-df_seqv[df_seqv[,i]==df_seqv[,number],]
    df_structure$aminoasid_number[i]<-nrow(v_leng)
  }
  protein<-colnames(df_seqv)[number]
  df_structure<-df_structure%>%filter(structure!=protein)
  df_structure<-df_structure%>%filter(aminoasid_number==max(aminoasid_number))
  write.csv(df_structure,"df_structure.csv",row.names = F)
  df_structure<-df_structure%>%mutate(pdb_name=NA)
  for (i in 1:nrow(df_structure)) {
    df_structure$pdb_name[i]<-strsplit(df_structure$structure[i],split = ":")[[1]][2]
  }
  df_structure<-df_structure%>%mutate(pdb=NA)
  for (i in 1:nrow(df_structure)) {
    df_structure$pdb[i]<-strsplit(df_structure$pdb_name[i],split = "_")[[1]][1]
  }

  for (i in 1:nrow(df_structure)) {
    print(paste0(main,"/structure_prediction/pdb/",df_structure$pdb[i],".pdb"))
    pdb<-read.pdb(paste0(main,"/structure_prediction/pdb/",df_structure$pdb[i],".pdb"))
    write.pdb(pdb,paste0(main,"/structure_prediction/pdb_fin/",df_structure$pdb[i],".pdb"))
  }
}
