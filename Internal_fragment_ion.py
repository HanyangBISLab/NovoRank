from numba import jit
from tqdm import tqdm
import pandas as pd
import numpy as np
import warnings
import os

warnings.filterwarnings(action='ignore')

amino_acid = {'G':57.02146372376, 'A':71.03711378804, 'S':87.03202841014, 'P':97.05276385232, 'V':99.0684139166,
              'T':101.04767847442, 'C':160.03064447804, 'L':113.08406398088, 'I':113.08406398088, 'N':114.04292744752,
              'D':115.02694303224, 'Q':128.0585775118, 'K':128.09496301826, 'E':129.04259309652, 'M':131.0404846066,
              'H':137.0589118628, 'F':147.0684139166, 'R':156.10111102874, 'Y':163.0633285387, 'W':186.07931295398,
              'm': 147.0353846066}

h20 = 18.01528
nh3 = 17.02655
proton = 1.00727647
neutron = 1.0086710869

fragment_tolerance = 0.025

def residue(sequence):

    b = []
    mass, mass_1, mass_2 = 0, 0, 0
    for j in sequence:
        mass += amino_acid[j]
        mass_1 = mass + neutron # isotope
        mass_2 = mass_1 + neutron # isotope
        b.append(mass)
        b.append(mass_1)
        b.append(mass_2)

    y = []
    mass, mass_1, mass_2 = 0, 0, 0
    for j in sequence[::-1]:
        mass += amino_acid[j]
        mass_1 = mass + neutron # isotope
        mass_2 = mass_1 + neutron # isotope
        y.append(mass)
        y.append(mass_1)
        y.append(mass_2)

    b = np.array(sorted(b))
    y = np.array(sorted(y))

    return b, y

@jit
def res_ion(b, y, charge, precur):

    temp_b = np.array([])
    for i in range(1, charge):
        temp_b = np.concatenate([temp_b, ((b+(proton*i))/i)])
        temp_b = np.concatenate([temp_b, ((b-h20+(proton*i))/i)])
        temp_b = np.concatenate([temp_b, ((b-nh3+(proton*i))/i)])

    temp_b = np.unique(temp_b)

    temp_y = np.array([])
    for i in range(1, charge):
        temp_y = np.concatenate([temp_y, ((y+h20+(proton*i))/i)])
        temp_y = np.concatenate([temp_y, ((y+h20-h20+(proton*i))/i)])
        temp_y = np.concatenate([temp_y, ((y+h20-nh3+(proton*i))/i)])

    temp_y = np.append(temp_y, np.array([precur]))
    temp_y = np.unique(temp_y)

    result = np.unique(np.concatenate([temp_b, temp_y]))

    return result

@jit
def remove_spec(spec_list, by):

    cnt = 0
    arr = []

    for i in by:

        if cnt == len(spec_list):
                break

        while True:
            if i-fragment_tolerance <= spec_list[cnt] and i+fragment_tolerance >= spec_list[cnt]:
                    arr.append(int(cnt))
                    cnt += 1
                    break

            if i+fragment_tolerance < spec_list[cnt]:
                break

            cnt += 1

            if cnt == len(spec_list):
                break

    temp = np.delete(spec_list, arr)

    return temp

def frag_ion(pep):

    pep_= pep[1:len(pep)-1]

    temp = []
    for n in range(1, len(pep_)+1):
        temp += [pep_[i:i+n]for i in range(len(pep_)-n+1)]

    temp = set(temp)

    temp_mass = set()
    append = temp_mass.add

    for i in temp:
        mass = 0
        for j in i:
            mass += amino_acid[j]
        append(mass)

    new_mass = np.array(sorted(temp_mass))

    return new_mass

@jit
def consider_charge(new_mass, charge):

    temp = np.array([])

    for i in range(1, charge):
        temp = np.concatenate([temp, ((new_mass+(proton*i))/i)])

    temp = np.unique(temp)

    return temp

@jit
def find_spec(spec_list, theo):

    temp = 0
    cnt = 0

    for i in theo:

        if cnt == len(spec_list):
                break
    
        while True:
            if i-fragment_tolerance <= spec_list[cnt] and i+fragment_tolerance >= spec_list[cnt]:
                    temp += 1
                    cnt += 1
                    break

            if i+fragment_tolerance < spec_list[cnt] :
                break

            cnt += 1
            if cnt == len(spec_list):
                break

    return temp

def internal_fragment_ion_cal(dataset, remove_path):

    file_list = os.listdir(remove_path)
    file_list.sort()

    temp = 0
    data = None
    spec_key = None
    internal_frag_ion = []
    append = internal_frag_ion.append
    
    for x, y, z, m, l in tqdm(dataset[['Source File', 'Scan number', 'Peptide', 'z', 'm/z']].values):

        d_fn = x.split('.')[0]
        d_scan = int(y)
        key = (d_fn, d_scan)

        if data == None:
            fn = file_list[temp]
            data = open(remove_path+'/'+fn).readlines()
            ind = 0

        if key == spec_key:

            b_ion, y_ion = residue(z)
            ion_list = res_ion(b_ion, y_ion, int(m), float(l))
            spec_arr = remove_spec(arr_mz, ion_list)
            residue_mass = frag_ion(z)
            final_mass = consider_charge(residue_mass, int(m))
            count = find_spec(spec_arr, final_mass)
            append(count)

            continue

        while True:

            if ind == len(data):
                temp += 1
                fn = file_list[temp]
                data = open(remove_path+'/'+fn).readlines()
                ind = 0

            if data[ind] == 'BEGIN IONS\n':
                cnt = 0
            else:
                cnt += 1
                if cnt == 2:
                    spec_scan = int(data[ind].split('.')[1])
                    spec_fn = data[ind].split('.')[0].split('=')[1]
                    spec_key = (spec_fn, spec_scan)

                    if d_fn > spec_fn:
                        temp += 1
                        fn = file_list[temp]
                        data = open(remove_path+'/'+fn).readlines()
                        ind = 0

                        continue

                    if key == spec_key:
                        arr_mz = np.array([])

                if key == spec_key and cnt > 5:
                    if data[ind] == 'END IONS\n' or data[ind] == 'END IONS':

                        b_ion, y_ion = residue(z)
                        ion_list = res_ion(b_ion, y_ion, int(m), float(l))
                        spec_arr = remove_spec(arr_mz, ion_list)
                        residue_mass = frag_ion(z)
                        final_mass = consider_charge(residue_mass, int(m))
                        count = find_spec(spec_arr, final_mass)
                        append(count)
                        ind += 1
                        break
                    else:
                        mz = float(data[ind].split()[0])
                        arr_mz = np.append(arr_mz, np.array([mz]))
            ind += 1
    
    dataset['internal_ion'] = internal_frag_ion
    
    return dataset

def seq_len(dataset):
    
    seq_len = []
    append = seq_len.append

    for i in dataset[['Sequence']].values:
        append(len(i[0]))

    dataset['seq_len'] = seq_len
    
    return dataset

def ifi_feature(dataset, remove_path):

    dataset_ = dataset[dataset['Peptide'].isnull()]
    dataset = dataset[dataset['Peptide'].notnull()]
    
    dataset = internal_fragment_ion_cal(dataset, remove_path)
    dataset = seq_len(dataset)
    
    dataset['int_frag_ion'] = dataset['internal_ion']/dataset['seq_len']
    dataset['int_frag_ion'] = np.log(dataset['int_frag_ion']+1)
    
    dataset_['internal_ion'] = np.nan
    dataset_['seq_len'] = np.nan
    dataset_['int_frag_ion'] = np.nan
    
    dataset = pd.concat([dataset, dataset_]).sort_values(by=['Source File', 'Scan number', 'Rank']).reset_index(drop=True)
    
    return dataset