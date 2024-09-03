import os
import pandas as pd
from tqdm import tqdm

def Spec_count(PATH):

    mgf_list = os.listdir(PATH)

    cnt = 0

    for i in mgf_list:
        data = open(PATH+'\\'+i).readlines()

        for j in data:
            if j == 'BEGIN IONS\n':
                cnt += 1

    print('Total Spectra : ', cnt)


def Extract_mgf_info(PATH):
    
    mgf_list = os.listdir(PATH)
    
    dataset = []
    append = dataset.append

    print("(MGF info) File index : ")
    
    for fn in tqdm(mgf_list):
        data = open(PATH+'\\'+fn).readlines()
    
        for i in data:
            if i == 'END IONS\n' or i == 'END IONS':
                d = fn, scan, charge, mass, rt
                append(d)
            else:
                if 'TITLE' in i:
                    scan = int(i.split('\n')[0].split('.')[1])
                elif 'CHARGE' in i:
                    charge = int(i.split('+\n')[0].split('=')[1])
                elif 'RTINSECONDS' in i:
                    rt = float(i.split('\n')[0].split('=')[1])
                elif 'PEPMASS' in i:
                    mass = float(i.split(' ')[0].split('=')[1])

    col = ['Source File', 'Scan number', 'z', 'm/z', 'rt']
    mgf_info = pd.DataFrame(dataset, columns=col)
    
    mgf_info['RT'] = mgf_info['rt']/60

    return mgf_info
