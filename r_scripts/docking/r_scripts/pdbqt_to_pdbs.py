# Using the magic encoding
# -*- coding: utf-8 -*-
import glob
import sys

def get_one_model (lines):
    res, res0 = [], []
    for ll in lines:
        if "MODEL" in ll:
            if len(res0)>0: res.append(res0)
            else: pass
            res0 = []
        if "ATOM" in ll: res0.append(ll)
        if "HETATM" in ll: res0.append(ll)
    else:
        if len(res0)>0: res.append(res0)
    
    return (res)
  
def work_onef (f_inp, save_name ):
    wf = open(f_inp, "r")
    lines = wf.readlines()
    wf.close()
    res = get_one_model (lines)
    ind = 1
    for model in res:
        save_name_ind = save_name % (ind)
        g = open(save_name_ind, "w")
        g.writelines(model)
        g.close()
        ind +=1

if __name__ == "__main__":
    args = sys.argv[1:]
    f_inp=args[0]
    save_name=args[1]
    work_onef (f_inp,save_name)
    

        
    


