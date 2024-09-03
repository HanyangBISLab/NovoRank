from DL_train import *
import os

def filtered_len_data(dataset):

    temp_x, temp_y = [], []

    for i, j in dataset[['Peptide_x','Peptide_y']].values:
        temp_x.append(len(i))
        temp_y.append(len(j))

    dataset['len_x'] = temp_x
    dataset['len_y'] = temp_y

    a = dataset['len_x'] <= 40
    b = dataset['len_y'] <= 40

    c = dataset['len_x'] > 40
    d = dataset['len_y'] > 40

    dataset_1 = dataset[a & b].reset_index(drop=True)
    dataset_2 = dataset[c | d].reset_index(drop=True)

    return dataset_1, dataset_2

def generator_test(new_data, path_dir):

    temp = 0  # file number
    length = 50000
    resolution = 0.1
    data = None
    spec_key = None

    file_list = os.listdir(path_dir)

    for fn, scan, seq_1, seq_2, z, score_x, score_y, delta_y, count, ion_x, ion_y, rt_x, rt_y, xcorr_x, xcorr_y, delta_ in \
            new_data[['Source File', 'Scan number', 'Peptide_x', 'Peptide_y', 'z_x', 'Score_x', 'Score_y', 'delta_y', 'new_count_x',
                      'int_frag_ion_x', 'int_frag_ion_y', 'diff_RT_x', 'diff_RT_y', 'xcorr_x', 'xcorr_y', 'delta_xcorr']].values:

        key = (fn.split('.')[0], int(scan))

        # data load
        if data == None:
            fn = file_list[temp]
            data = open(path_dir + '/' + fn).readlines()
            ind = 0

        while True:

            if ind == len(data):
                temp += 1

                fn = file_list[temp]
                data = open(path_dir + '/' + fn).readlines()
                ind = 0

            if data[ind] == 'BEGIN IONS\n':
                cnt = 0
            else:
                cnt += 1
                if cnt == 2:
                    spec_scan = int(data[ind].split('.')[1])
                    spec_fn = data[ind].split('.')[0].split('=')[1]
                    spec_key = (spec_fn, spec_scan)

                    if fn.split('.')[0] > spec_fn:
                        temp += 1
                        fn = file_list[temp]
                        data = open(path_dir + '/' + fn).readlines()
                        ind = 0

                        continue

                if key == spec_key and cnt == 4:
                    spectrum = make_spectrum(length).astype('float64')

                if key == spec_key and cnt > 5:
                    if data[ind] == 'END IONS\n' or data[ind] == 'END IONS':

                        spectrum = normalization(spectrum)
                        seq_1 = sequence_encoding(seq_1, z)
                        seq_2 = sequence_encoding(seq_2, z)
                        add_1 = add_psm(score_x, ion_x, rt_x, xcorr_x)
                        add_2 = add_psm(score_y, ion_y, rt_y, xcorr_y)
                        add_add = add_info(delta_y, count, delta_)

                        yield {'input_spec': spectrum, 'input_seq_1': seq_1, 'input_seq_2': seq_2,
                               'input_add_1': add_1,
                               'input_add_2': add_2, 'input_add_add': add_add}

                        ind += 1
                        break
                    else:
                        mz = float(data[ind].split()[0])
                        intensity = float(data[ind].split()[1].split('\n')[0])
                        loc = int(mz // resolution)

                        if loc > 50000:
                            pass
                        else:
                            spectrum[loc] = sum_intensity(spectrum[loc], intensity)

            ind += 1

def pred_zero_to_one(dataset):

    temp = []
    append = temp.append

    for i in dataset[['predict']].values:
        if i[0] <= 0.5:
            append(0)
        else:
            append(1)

    return temp

def test_results_dataset(dataset):

    temp = []
    append = temp.append

    for i, pep_1, pep_2 in dataset[['pred', 'Peptide_x', 'Peptide_y']].values:
        if i == 0:
            append(pep_1)
        else:
            append(pep_2)

    return temp