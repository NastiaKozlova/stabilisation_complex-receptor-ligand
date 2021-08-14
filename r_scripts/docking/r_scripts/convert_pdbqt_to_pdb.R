part_start <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
setwd(part_start)
if (!dir.exists(paste0(part_start,"analysis/"))){dir.create(paste0(part_start,"analysis/"))}
df_convert<-paste0('# Using the magic encoding\n',
                   '# -*- coding: utf-8 -*-\n',
                   'import glob\n\n',
                   'def get_one_model (lines):\n',
                   '    res, res0 = [], []\n',
                   '    for ll in lines:\n',
                   '        if "MODEL" in ll:\n',
                   '            if len(res0)>0: res.append(res0)\n',
                   '            else: pass\n',
                   '            res0 = []\n',
                   '        if "ATOM" in ll: res0.append(ll)\n',
                   '    else:\n',
                   '        if len(res0)>0: res.append(res0)\n',
                   '    return (res)\n\n',
                   'def gey_name(f, f_save, ind):\n',
                   '    li = f.split(r"/") # В Windows другой слеш!!!\n',
                   '    way = r"/".join(li[:-1]) \n',
                   '    name = ".".join(li[-1].split(".")[:-1])\n',
                   '    save_name = "%s/%s_MODEL_%d.pdb" % (f_save, name, ind)\n',
                   '    return (save_name)\n',
                   'def open_read (f_open, f_save):\n',
                   '    ff = glob.glob(f_open)\n',
                   '    for f in ff:\n',
                   '        wf = open(f, "r")\n',
                   '        lines = wf.readlines()\n',
                   '        wf.close()\n',
                   '        res = get_one_model (lines)\n',
                   '        ind = 1\n',
                   '        for model in res:\n',
                   '            save_name = gey_name(f,f_save, ind)\n',
                   '            g = open(save_name, "w")\n',
                   '            g.writelines(model)\n',
                   '            g.close()\n',
                   '            ind +=1\n\n',
                   'f_open = "',part_start,'out/*.pdbqt"\n',
                   'f_save = "',part_start,'analysis/"\n',
                   'open_read (f_open, f_save)\n')
write.table(df_convert,paste0(part_start,"convert_pdbqt_to_pdb.py"),quote = F,na = "",col.names = F,row.names = F)


