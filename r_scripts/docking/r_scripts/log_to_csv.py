import glob

folder_log = 'ACHE/docking/docking_first/log/*.log'
folder_svg = 'ACHE/docking/docking_first/din/log/%s.csv'

def convert_one_file (fname, fsave):
    f = open(fname, "r")
    text = f.read()
    f.close()
    text0 = text.split("mode |")
    w_rows = []
    if len(text0)>0:
        rows = text0[1].split("Writing")
        if len(rows)>0:
            need_rows = rows[0].split("\n")
            w_rows = [",".join(x.split()) + "\n" for x in need_rows[3:]]
    if len(w_rows) > 0:
        g = open (fsave, "w")
        g.writelines (w_rows)
        g.close()
        print ("file %s writen" % (fsave))
    else:
        print ("file %s NO_writen" % (fsave))


def convert_logs_to_svg (folder_log, folder_svg):
    fnames = glob.glob(folder_log)
    for fname in fnames:
        fsave = folder_svg % (fname.split("/")[-1][:-4])
        convert_one_file (fname, fsave)


convert_logs_to_svg (folder_log, folder_svg)
