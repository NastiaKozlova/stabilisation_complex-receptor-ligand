library(bio3d)
library(dplyr)
part_start = commandArgs(trailingOnly=TRUE)
setwd(part_start)
if (!dir.exists(paste0("pdb_fin"))){dir.create(paste0("pdb_fin"))}
seq<-read.fasta(paste0(part_start,"protein_sequence.fasta"))
fasta<-read.fasta(paste0(part_start,"seqs.fasta"))
tost<-seqbind(seq,fasta)
aln<-seqaln(tost,seqgroup = T, exefile="muscle") 
df_seq<-as.data.frame(t((aln$ali)),stringsAsFactors = F)
v_names<-colnames(df_seq)
df_structure<-data.frame(matrix(ncol=2,nrow = length(v_names)))
colnames(df_structure)<-c("structure","aminoasid_number")
df_structure$structure<-v_names
df_seq[df_seq=="-"]<-NA
i<-1
for (i in 1:nrow(df_structure)) {
  if(df_structure$structure[i]!=seq$id){
    df_seq_TEMP<-df_seq%>%select(c(df_structure$structure[i],seq$id))
    colnames(df_seq_TEMP)<-c("tested_seq","target_seq")
    df_seq_TEMP<-df_seq_TEMP%>%filter(!is.na(target_seq))
    df_seq_TEMP<-df_seq_TEMP%>%filter(tested_seq==target_seq)
    df_structure$aminoasid_number[i]<-nrow(df_seq_TEMP)
  }
}
df_structure<-df_structure%>%filter(!is.na(aminoasid_number))
df_structure<-df_structure%>%filter(aminoasid_number==max(aminoasid_number))
write.csv(df_structure,"df_structure.csv",row.names = F)
df_structure<-df_structure%>%mutate(pdb_name=NA)
for (i in 1:nrow(df_structure)) {
  df_structure$pdb_name[i]<-strsplit(df_structure$structure[i],split = "_")[[1]][1]
}
i<-1
for (i in 1:nrow(df_structure)) {
  a<-strsplit(df_structure$pdb_name[i],split = ":",fixed = T)[[1]]
  if(a[1]=="PDB"){
    df_structure$pdb_name[i]<-a[2]
  }
}
for (i in 1:nrow(df_structure)) {
  print(paste0("pdb/",df_structure$pdb_name[i],".pdb"))
  if(file.exists(paste0("pdb/", df_structure$pdb_name[i], ".pdb"))){
    pdb<-read.pdb(paste0("pdb/", df_structure$pdb_name[i], ".pdb"))
    write.pdb(pdb,paste0("pdb_fin/", df_structure$pdb_name[i], ".pdb"))
  }
}

