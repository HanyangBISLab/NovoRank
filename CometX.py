import os
from tqdm import tqdm

def CometX_mgf(dataset, mgf_dir, save_path):

    dataset_new = dataset[dataset['Peptide'].notnull()]

    mgf_list = os.listdir(mgf_dir)

    scan = -1
    file_index = 0

    for ind, (i, j, k) in enumerate(tqdm(dataset_new[['Source File', 'Scan number', 'Peptide']].values)):

        if ind == 0:
            spec = open(mgf_dir+'\\'+mgf_list[file_index]).readlines()
            loc = 0
            f = open(save_path+'\\'+mgf_list[file_index], 'w')

        if scan == j:
            for idx, row in enumerate(temp):
                if idx == 2:
                    seq = k.replace('m', 'M+15.9949').replace('C', 'C+57.021464')
                    f.write('SEQ='+seq+'\n')
                f.write(row)
        else:
            if loc == len(spec):
                f.close()
                file_index += 1
                if file_index > len(mgf_list):
                    break
                spec = open(mgf_dir+'\\'+mgf_list[file_index]).readlines()
                loc = 0
                f = open(save_path+'\\'+mgf_list[file_index], 'w')

            temp = []
            append = temp.append

            while True:

                if loc == len(spec):
                    f.close()
                    file_index += 1
                    if file_index > len(mgf_list):
                        break
                    spec = open(mgf_dir+'\\'+mgf_list[file_index]).readlines()
                    loc = 0
                    f = open(save_path+'\\'+mgf_list[file_index], 'w')

                if spec[loc] == 'BEGIN IONS\n':
                    start = loc
                    append(spec[loc])
                elif spec[loc] == 'END IONS\n' or spec[loc] == 'END IONS':
                    append(spec[loc])

                    if scan == j:
                        for idx, row in enumerate(temp):
                            if idx == 2:
                                seq = k.replace('m', 'M+15.9949').replace('C', 'C+57.021464')
                                f.write('SEQ='+seq+'\n')
                            f.write(row)
                        loc += 1
                        break
                    else:
                        temp = []
                        append = temp.append
                elif loc == start + 1:
                    append(spec[loc])
                    scan = int(spec[loc].split('=')[1].split('\n')[0])
                else:
                    append(spec[loc])
                loc += 1

    f.close()