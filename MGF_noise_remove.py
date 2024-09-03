import os
from tqdm import tqdm

def Remove_noise(temp_list):
    
    temp = []
    for i in temp_list:
        temp.append(float(i.split()[1].split('\n')[0]))
    
    while True:
        idx = temp.index(min(temp))
        temp_list.remove(temp_list[idx])
        temp.remove(temp[idx])

        if len(temp_list) < 11:
            break
        
    return temp_list

def Make_remove_mgf(path_dir, save_dir):

    file_list = os.listdir(path_dir)

    print('(Remove Noise) File index : ')
    for i in tqdm(file_list):
        with open(path_dir+'/'+i, 'r') as infile:
            data = infile.readlines()

        with open(save_dir+"\\"+i.split('.')[0]+'_remove'+'.mgf', 'w') as outfile:
            ind = 0
            for i in data:
                if i == 'BEGIN IONS\n':
                    outfile.write(i)
                    ind += 1
                    cnt = 0
                    temp = []
                    window_min = 0
                    window_max = 100

                elif i == 'END IONS\n' or i == 'END IONS':
                    if len(temp) > 10:
                        temp = Remove_noise(temp)
                        for j in temp:
                            outfile.write(j)
                    else:
                        for j in temp:
                            outfile.write(j)
                    outfile.write(i)

                else:
                    cnt += 1
                    if cnt > 5:
                        if float(i.split()[0]) > window_max:
                            while float(i.split()[0]) > window_max:
                                window_min += 100
                                window_max += 100
                            if len(temp) == 0:
                                temp.append(i)
                                continue

                            elif len(temp) > 10:
                                temp = Remove_noise(temp)
                                for j in temp:
                                    outfile.write(j)
                                temp = []

                            else:
                                for j in temp:
                                    outfile.write(j)
                                temp = []
                            temp.append(i)

                        else:
                            temp.append(i)

                    else:
                        outfile.write(i)