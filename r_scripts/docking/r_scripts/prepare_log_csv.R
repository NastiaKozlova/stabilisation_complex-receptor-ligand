part_name <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
setwd(part_name)
if (!dir.exists(paste0(part_name,"din/"))){dir.create(paste0(part_name,"din/"))}
if (!dir.exists(paste0(part_name,"din/log/"))){dir.create(paste0(part_name,"din/log/"))}
df_convert<-paste0('import glob\n\n',
                   'folder_log = "',part_name,'log/*.log"\n',
                   'folder_svg = "',part_name,'din/log/%s.csv"\n\n',
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
                   'def convert_logs_to_svg (folder_log, folder_svg):\n',
                   '    fnames = glob.glob(folder_log)\n',
                   '    for fname in fnames:\n',
                   '        fsave = folder_svg % (fname.split("/")[-1][:-4])\n',
                   '        convert_one_file (fname, fsave)\n',
                   'convert_logs_to_svg (folder_log, folder_svg)')
write.table(df_convert,paste0(part_name,"prepare_log_csv.py"),quote = F,na = "",col.names = F,row.names = F)
