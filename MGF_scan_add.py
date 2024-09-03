import os
from tqdm import tqdm

def Split_scan(mgf):
    
    temp = []
    append = temp.append
    
    for idx, i in enumerate(mgf):
        if i == 'BEGIN IONS\n':
            start = idx
        elif 'TITLE' in i:
            scan = int(i.split('\n')[0].split('.')[1])
        else:
            if i == 'END IONS\n' or i == 'END IONS':
                end = idx + 1
                s = []
                spec = mgf[start:end]
                for ind, i in enumerate(spec):
                    if ind == 0:
                        s.append(i)
                        s.append('SCANS='+str(scan)+"\n")
                    else:
                        s.append(i)
                append(tuple(s))
    
    return temp

def Write(MGF_PATH, file_name, mgf):
    
    f = open(MGF_PATH+"\\"+file_name, 'w')
    
    for i in mgf:
        for j in i:
            f.write(j)
    f.close()

def Mgf(path_dir):
        
    mgf_list = os.listdir(path_dir)

    print("(ADD SCAN) File index : ")
    for fn in tqdm(mgf_list):
        contents = open(path_dir+'\\'+fn).readlines()

        if 'SCANS=' in contents[1]:
            pass
        else:
            add_scan = Split_scan(contents)
            Write(path_dir, fn, add_scan)