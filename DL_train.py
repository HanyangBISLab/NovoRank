from tensorflow.keras.layers import Conv1D, LSTM, Bidirectional, Dense, Flatten # BatchNormalization
from tensorflow.keras.layers import MaxPool1D, LeakyReLU, Dropout, Reshape, concatenate, Activation
from tensorflow import keras
import tensorflow as tf

from sklearn.model_selection import GroupShuffleSplit
from numba import jit
import numpy as np
import os

@jit
def make_aa(element):
    return np.zeros(element, dtype=np.int8)


@jit
def sequence_encoding(sequence, z):
    charge = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5}

    aa = {'A': 6, 'C': 7, 'D': 8, 'E': 9, 'F': 10, 'G': 11, 'H': 12, 'I': 13, 'K': 14, 'L': 15, 'M': 16, 'm': 16,
          'N': 17, 'P': 18, 'Q': 19, 'R': 20, 'S': 21, 'T': 22, 'V': 23, 'W': 24, 'Y': 25}

    modification = {'m': 26, 'C': 27}

    temp = []  # peptide
    modi = False  # init

    for i in sequence:
        temp_np = make_aa(28)
        temp_np[aa[i]] = 1
        temp_np[charge[int(z)]] = 1

        if i == 'C':
            temp_np[modification[i]] = 1
        elif i == 'm':
            temp_np[modification[i]] = 1

        temp.append(temp_np)

    while len(temp) < 40:
        temp_np = make_aa(28)
        temp.append(temp_np)

    arr = np.array(temp)

    return arr


@jit
def normalization(spectrum):
    factor = max(spectrum)
    spectrum = spectrum / factor

    return spectrum


@jit
def lebel_append(label):
    new_label = np.array([label])

    return new_label

@jit
def add_psm(score, ion, rt, xcorr):
    new_info = np.array([score, ion, rt, xcorr])

    return new_info

@jit
def add_info(delta, count, delta_):
    new_info = np.array([delta, count, delta_])

    return new_info

@jit
def make_spectrum(length):
    spec = np.zeros((length))

    return spec

@jit
def sum_intensity(cur, new):
    return cur + np.array(new)


def used_dataset(dataset):
    temp = []
    append = temp.append

    for idx, (i, j) in enumerate(dataset[['Label_x', 'Label_y']].values):
        if i + j == 1:
            append(idx)

    dataset = dataset.loc[temp]

    return dataset

def train_val_split(dataset, val_size):
    groups = dataset['GT_x']
    gss = GroupShuffleSplit(n_splits=1, test_size=val_size)

    for train_idx, val_idx in gss.split(dataset, groups=groups):
        print("TRAIN_DATA_LENGTH :", len(train_idx))
        print("VAL_DATA_LENGTH :", len(val_idx))

    train_data = dataset.iloc[train_idx].reset_index(drop=True)
    val_data = dataset.iloc[val_idx].reset_index(drop=True)

    return train_data, val_data

def generator_train(train_data, path_dir):

    temp = 0  # file number
    length = 50000
    resolution = 0.1
    data = None
    spec_key = None

    file_list = os.listdir(path_dir)

    for fn, scan, seq_1, seq_2, z, score_x, score_y, delta_y, count, ion_x, ion_y, rt_x, rt_y, xcorr_x, xcorr_y, delta_, label in \
        train_data[['Source File', 'Scan number', 'Peptide_x', 'Peptide_y', 'z_x', 'Score_x', 'Score_y', 'delta_y', 'new_count_x',
                    'int_frag_ion_x', 'int_frag_ion_y', 'diff_RT_x', 'diff_RT_y', 'xcorr_x', 'xcorr_y', 'delta_xcorr',
                    'Label_x']].values:

        key = (fn.split('.')[0], int(scan))

        # data load
        if data == None:
            fn = file_list[temp]
            data = open(path_dir + '/' + fn).readlines()
            ind = 0

        #         if key == spec_key:

        #             spectrum = normalization(spectrum)
        #             seq_1 = sequence_encoding(seq_1, z)
        #             seq_2 = sequence_encoding(seq_2, z)
        #             add_1 = add_info(score_x, delta_x, count, ion_x, rt_x)
        #             add_2 = add_info(score_y, delta_y, count, ion_y, rt_y)
        #             label_ = lebel_append(label)

        #             yield  {'input_spec': spectrum, 'input_seq_1': seq_1 , 'input_seq_2': seq_2 , 'input_add_1': add_1, 'input_add_2': add_2}, label_

        #             continue

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
                    spectrum = make_spectrum(length)

                if key == spec_key and cnt > 5:
                    if data[ind] == 'END IONS\n' or data[ind] == 'END IONS':

                        spectrum = normalization(spectrum)
                        seq_1 = sequence_encoding(seq_1, z)
                        seq_2 = sequence_encoding(seq_2, z)
                        add_1 = add_psm(score_x, ion_x, rt_x, xcorr_x)
                        add_2 = add_psm(score_y, ion_y, rt_y, xcorr_y)
                        add_add = add_info(delta_y, count, delta_)
                        label_ = lebel_append(label)

                        yield {'input_spec': spectrum, 'input_seq_1': seq_1, 'input_seq_2': seq_2, 'input_add_1': add_1,
                               'input_add_2': add_2, 'input_add_add': add_add}, label_

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


def generator_val(val_data, path_dir):

    temp = 0  # file number
    length = 50000
    resolution = 0.1
    data = None
    spec_key = None

    file_list = os.listdir(path_dir)

    for fn, scan, seq_1, seq_2, z, score_x, score_y, delta_y, count, ion_x, ion_y, rt_x, rt_y, xcorr_x, xcorr_y, delta_, label in \
        val_data[['Source File', 'Scan number', 'Peptide_x', 'Peptide_y', 'z_x', 'Score_x', 'Score_y', 'delta_y', 'new_count_x',
                  'int_frag_ion_x', 'int_frag_ion_y', 'diff_RT_x', 'diff_RT_y', 'xcorr_x', 'xcorr_y', 'delta_xcorr',
                  'Label_x']].values:

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
                    spectrum = make_spectrum(length)

                if key == spec_key and cnt > 5:
                    if data[ind] == 'END IONS\n' or data[ind] == 'END IONS':

                        spectrum = normalization(spectrum)
                        seq_1 = sequence_encoding(seq_1, z)
                        seq_2 = sequence_encoding(seq_2, z)
                        add_1 = add_psm(score_x, ion_x, rt_x, xcorr_x)
                        add_2 = add_psm(score_y, ion_y, rt_y, xcorr_y)
                        add_add = add_info(delta_y, count, delta_)
                        label_ = lebel_append(label)

                        yield {'input_spec': spectrum, 'input_seq_1': seq_1, 'input_seq_2': seq_2, 'input_add_1': add_1,
                               'input_add_2': add_2, 'input_add_add': add_add}, label_

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


def NovoRank_model():
    spec = keras.layers.Input(shape=(50000, 1))
    seq = keras.layers.Input(shape=(40, 28))
    plus = keras.layers.Input(shape=(4))

    encoder = Conv1D(8, 30, strides=1, padding="valid")(spec)
    encoder = LeakyReLU(alpha=0.01)(encoder)
    encoder = Dropout(0.1)(encoder)
    encoder = MaxPool1D(pool_size=30, strides=30, padding='valid')(encoder)
    encoder = Conv1D(16, 30, strides=1, padding="valid")(encoder)
    encoder = LeakyReLU(alpha=0.01)(encoder)
    encoder = Dropout(0.1)(encoder)
    encoder = MaxPool1D(pool_size=30, strides=30, padding='valid')(encoder)
    encoder = Flatten()(encoder)
    encoder = Dense(16)(encoder)
    encoder = LeakyReLU(alpha=0.01)(encoder)
    encoder = Dropout(0.1)(encoder)
    encoder = Dense(16)(encoder)
    encoder = LeakyReLU(alpha=0.01)(encoder)
    encoder = Dropout(0.1)(encoder)
    encoder = Dense(16)(encoder)

    bilstm = Bidirectional(LSTM(8), merge_mode='concat')(seq)
    bilstm = Dense(16)(bilstm)
    bilstm = LeakyReLU(alpha=0.01)(bilstm)
    bilstm = Dropout(0.1)(bilstm)
    bilstm = Dense(16)(bilstm)
    bilstm = LeakyReLU(alpha=0.01)(bilstm)
    bilstm = Dropout(0.1)(bilstm)
    bilstm = Dense(16)(bilstm)

    output = concatenate([encoder, bilstm, plus])
    output = Dense(32)(output)
    output = LeakyReLU(alpha=0.01)(output)
    output = Dropout(0.1)(output)
    output = Dense(32)(output)
    output = LeakyReLU(alpha=0.01)(output)
    output = Dropout(0.1)(output)
    output = Dense(16)(output)

    output = Reshape((-1,))(output)

    model_temp = tf.keras.Model(inputs=[spec, seq, plus], outputs=output)

    return model_temp

def psm_model(model_1, model_2, concat_add):

    psm = concatenate([model_1, model_2, concat_add])

    output_ = Dense(32)(psm)
    output_ = LeakyReLU(alpha=0.01)(output_)
    output_ = Dropout(0.1)(output_)
    output_ = Dense(32)(output_)
    output_ = LeakyReLU(alpha=0.01)(output_)
    output_ = Dropout(0.1)(output_)
    output_ = Dense(1)(output_)
    output_ = Activation('sigmoid')(output_)

    return output_