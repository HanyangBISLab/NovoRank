import os
import pandas as pd
from tqdm import tqdm

def cross_corrlation_info(path_dir):

    fileEx = '.tsv'
    file_list = [file for file in os.listdir(path_dir) if file.endswith((fileEx))]

    temp = []
    append = temp.append

    for j in tqdm(file_list):

        spec = open(path_dir+'\\'+j).readlines()

        for idx, i in enumerate(spec):
            if idx == 0:
                pass
            else:
                fn = i.split('.')[0]+'.mgf'
                scan = i.split('.')[1]
                charge = i.split('\t')[1].split(' ')[0]
                seq = i.split('\t')[4].replace('M+15.9949', 'm').replace('C+57.021464', 'C')
                xcorr = i.split('\t')[5].split('\n')[0]
                data = (fn, scan, int(charge), seq, float(xcorr))
                append(data)

    columns = ['Source File', 'Scan number', 'z', 'Peptide', 'xcorr']
    df = pd.DataFrame(temp, columns=columns).astype({'Scan number': 'int'})

    return df

def remove_duplication(dataset):
    dataset = dataset.sort_values(by=['Source File', 'Scan number', 'Rank'])

    dataset1 = dataset.drop_duplicates(subset=['Source File', 'Scan number'], keep='first')
    dataset2 = dataset.drop_duplicates(subset=['Source File', 'Scan number'], keep='last')

    new = pd.merge(dataset1, dataset2, how='outer', on=['Source File', 'Scan number']).reset_index(drop=True)

    return new