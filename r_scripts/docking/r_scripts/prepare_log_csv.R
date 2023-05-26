part_analysis <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
setwd(part_analysis)
if (!dir.exists(paste0(part_analysis,"din/"))){dir.create(paste0(part_analysis,"din/"))}
if (!dir.exists(paste0(part_analysis,"din/convertion_python"))){dir.create(paste0(part_analysis,"din/convertion_python"))}
if (!dir.exists(paste0(part_analysis,"din/log/"))){dir.create(paste0(part_analysis,"din/log/"))}
log_files<-list.files(paste0(part_analysis,"log/"))
df_convert<-data.frame(matrix(nrow = length(log_files),ncol = 5))
colnames(df_convert)<-c("start_name","save_name","converted","script_name","file_size")
df_convert$start_name<-log_files
i<-1
for (i in 1:nrow(df_convert)) {
   
   a<-strsplit(df_convert$start_name[i],split = ".",fixed = T)[[1]][1]
   df_convert$script_name[i]<-paste0(a)
   df_convert$save_name[i]<-paste0(part_analysis,"din/log/",a,".csv")
   df_convert$file_size[i]<-file.size(paste0(part_analysis,"log/",df_convert$start_name[i]))
   df_convert$converted[i]<-file.exists(paste0(part_analysis,"din/log/",a,".csv"))
}
df_convert<-df_convert%>%mutate(start_name=paste0(part_analysis,"log/",script_name,".log"))
df_convert<-df_convert%>%filter(!converted)
df_convert<-df_convert%>%filter(file_size>1200)

i<-1
if (nrow(df_convert)>0) {
   
   for (i in 1:nrow(df_convert)) {
      
      
      df_script<-paste0('import glob\n\n',
                        'start_name = "',df_convert$start_name[i],'"\n',
                        'save_name = "',df_convert$save_name[i],'"\n\n',
                        'def convert_one_file (fname, fsave):\n',
                        '    f = open(fname, "r")\n',
                        '    text = f.read()\n',
                        '    f.close()\n',
                        '    text0 = text.split("mode |")\n',
                        '    w_rows = []\n',
                        '    if len(text0)>0:\n',
                        '        rows = text0[1].split("Writing")\n',
                        '        if len(rows)>0:\n',
                        '            need_rows = rows[0].split("\\n")\n',
                        '            w_rows = [",".join(x.split()) + "\\n" for x in need_rows[3:]]\n',
                        '    if len(w_rows) > 0:\n',
                        '        g = open (fsave, "w")\n',
                        '        g.writelines (w_rows)\n',
                        '        g.close()\n',
                        '        print ("file %s writen" % (fsave))\n',
                        '    else:\n',
                        '        print ("file %s NO_writen" % (fsave))\n',
                        'convert_one_file (start_name, save_name)')
      write.table(df_script,paste0(part_analysis,"din/convertion_python/",df_convert$script_name[i],".py"),quote = F,na = "",col.names = F,row.names = F)
      system(command = paste0("python3 ",part_analysis,"din/convertion_python/",df_convert$script_name[i],".py"),ignore.stdout=T,wait = T)
      
   }
}
#print(paste0(df_convert$start_name[i],".log"))
