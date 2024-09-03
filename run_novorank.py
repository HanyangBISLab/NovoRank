import os
import warnings
import argparse

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

warnings.filterwarnings(action='ignore')

from Config import *
from Xcorr import *
from NEW_candidates import *
from DL_test import *

parser = argparse.ArgumentParser(description='Description')
parser.add_argument('config', help='config text file for run_novorank.py')

args = parser.parse_args()

configs = read_config(args.config)

datasets = dataload(configs['RESULT_NAME'])

# After Xcorr calculation
cometX_results_path = configs['XCORR_RESULT']
xcorr_info = cross_corrlation_info(cometX_results_path)

new_datasets = pd.merge(datasets, xcorr_info, on=['Source File', 'Scan number', 'Peptide', 'z'], how='outer')
# new_datasets['Xcorr'] = np.log(1 + new_datasets['xcorr'])
new_datasets = new_datasets[new_datasets['Peptide'].notnull()].reset_index(drop=True)

remain_datasets = new_datasets[new_datasets['xcorr'].isnull()].reset_index(drop=True)
new_datasets = new_datasets[new_datasets['xcorr'].notnull()].reset_index(drop=True)

new_df = remove_duplication(new_datasets)
new_df['delta_xcorr'] = new_df['xcorr_x'] - new_df['xcorr_y']

# Deep Learning
new_data, remain_datasets_2 = filtered_len_data(new_df)
# print(tf.__version__)
# print(tf.test.is_built_with_cuda())
# # print(tf.sysconfig.get_build_info())
gpus = tf.config.experimental.list_physical_devices('GPU')
# print(gpus)
# tf.config.experimental.set_visible_devices(gpus[0], 'GPU')

if gpus:
    try:
        tf.config.experimental.set_memory_growth(gpus[0], True)
        os.environ["CUDA_VISIBLE_DEVICES"] = "0"
    except RuntimeError as e:
        print(e)
else:
    os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

if configs['TRAIN'] == True:
    new_data_ = used_dataset(new_data)
    train, val = train_val_split(new_data_, configs['VAL_SIZE'])

    train = train[['Source File', 'Scan number', 'Peptide_x', 'Peptide_y', 'z_x', 'Rank_x', 'Rank_y', 'Score_x', 'Score_y',
                   'delta_y', 'new_count_x', 'int_frag_ion_x', 'int_frag_ion_y', 'diff_RT_x', 'diff_RT_y', 'xcorr_x', 'xcorr_y',
                   'delta_xcorr', 'Label_x']].copy()
    val = val[['Source File', 'Scan number', 'Peptide_x', 'Peptide_y', 'z_x', 'Rank_x', 'Rank_y', 'Score_x', 'Score_y',
               'delta_y', 'new_count_x', 'int_frag_ion_x', 'int_frag_ion_y', 'diff_RT_x', 'diff_RT_y', 'xcorr_x', 'xcorr_y',
               'delta_xcorr', 'Label_x']].copy()

    output_shape = ({'input_spec': [50000], 'input_seq_1': [40, 28], 'input_seq_2': [40, 28],
                     'input_add_1': [4], 'input_add_2': [4], 'input_add_add': [3]}, [1])

    data_train = tf.data.Dataset.from_generator(lambda: generator_train(train, configs['MGF_PATH']),
                                                output_types=({'input_spec': tf.float64, 'input_seq_1': tf.float64, 'input_seq_2': tf.float64,
                                                               'input_add_1': tf.float64, 'input_add_2': tf.float64, 'input_add_add': tf.float64},
                                                              tf.float64),
                                                output_shapes=output_shape)
    data_val = tf.data.Dataset.from_generator(lambda: generator_val(val, configs['MGF_PATH']),
                                              output_types=({'input_spec': tf.float64, 'input_seq_1': tf.float64, 'input_seq_2': tf.float64,
                                                             'input_add_1': tf.float64, 'input_add_2': tf.float64, 'input_add_add': tf.float64},
                                                            tf.float64),
                                              output_shapes=output_shape)

    data_train = data_train.batch(configs['BATCH'])
    data_val = data_val.batch(configs['BATCH'])

    spectra = keras.layers.Input(shape=(50000, 1), name='input_spec')
    sequence_1 = keras.layers.Input(shape=(40, 28), name='input_seq_1')
    sequence_2 = keras.layers.Input(shape=(40, 28), name='input_seq_2')
    concat_1 = keras.layers.Input(shape=(4), name='input_add_1')
    concat_2 = keras.layers.Input(shape=(4), name='input_add_2')
    concat_add = keras.layers.Input(shape=(3), name='input_add_add')

    base_model = NovoRank_model()

    model_1 = base_model([spectra, sequence_1, concat_1])
    model_2 = base_model([spectra, sequence_2, concat_2])

    output = psm_model(model_1, model_2, concat_add)

    model = tf.keras.Model(inputs=[spectra, sequence_1, sequence_2, concat_1, concat_2, concat_add], outputs=output)
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics='accuracy')

    # model.summary()
    print('Train start')
    print()
    history = model.fit(data_train, validation_data=data_val, epochs=2, verbose=1) # configs['EPOCH']

    model.save(configs['MODEL_NAME'])
    print()
    print('Train done')
else:
    pre_trained_model = tf.keras.models.load_model(configs['PRE_TRAINED_MODEL'])

    output_shape = ({'input_spec': [50000], 'input_seq_1': [40, 28], 'input_seq_2': [40, 28],
                     'input_add_1': [4], 'input_add_2': [4], 'input_add_add': [3]})
    data_test = tf.data.Dataset.from_generator(lambda: generator_test(new_data, configs['MGF_PATH']),
                                               output_types=({'input_spec': tf.float64, 'input_seq_1': tf.float64, 'input_seq_2': tf.float64,
                                                              'input_add_1': tf.float64, 'input_add_2': tf.float64, 'input_add_add': tf.float64}),
                                               output_shapes=output_shape)

    data_test = data_test.batch(configs['BATCH'])

    print('Test start')
    print()
    output = pre_trained_model.predict(data_test, verbose=1)

    new_data['predict'] = output
    new_data['pred'] = pred_zero_to_one(new_data)

    new_data['Peptide'] = test_results_dataset(new_data)

    results_1 = remain_datasets[['Source File', 'Scan number', 'Peptide']]
    results_2 = remain_datasets_2[['Source File', 'Scan number', 'Peptide_x']].rename(columns={'Peptide_x': 'Peptide'})
    results_3 = new_data[['Source File', 'Scan number', 'Peptide']]

    results_final = pd.concat([results_1, results_2, results_3]).sort_values(by=['Source File', 'Scan number']).reset_index(drop=True)

    results_final.to_csv(configs['FINAL_RESULT'], index=False)
    print()
    print('Test done')