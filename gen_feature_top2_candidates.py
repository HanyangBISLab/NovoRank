import shutil
import argparse

from Config import *
from MGF_scan_add import *
from Cluster_csv import *
from MGF_info import *
from MGF_noise_remove import *
from NEW_candidates import *
from Internal_fragment_ion import *
from Cal_rt_feature import *
from CometX import *

# pd.options.display.max_columns = 100

# Press the green button in the gutter to run the script.

parser = argparse.ArgumentParser(description='Description')
parser.add_argument('config', help='config text file for gen_feature_top2_candidates.py')

if __name__ == '__main__':

    args = parser.parse_args()

    # read config file
    print('Read config file')
    config = read_config(args.config)
    print()

    # pre-processing
    print('< Pre-processing >')
    Mgf(config['MGF_PATH'])

    cluster_list = os.listdir(config['CLUSTER_PATH'])
    mgf_list = os.listdir(config['MGF_PATH'])
    index_dic = File_index(mgf_list)

    cluster_info = Cluster_csv(config['CLUSTER_PATH'], cluster_list, index_dic)

    Spec_count(config['MGF_PATH'])
    mgf_info = Extract_mgf_info(config['MGF_PATH'])

    remove_path = config['MGF_REMOVE']
    if os.path.exists(remove_path):
        shutil.rmtree(remove_path)
    os.mkdir(remove_path)
    Make_remove_mgf(config['MGF_PATH'], remove_path) # 노이즈 제거한 mgf : 100Da 상위 10개 픽

    # making candidates
    print('< Making candidates >')

    print('>> Dataload ...')
    de_novo = dataload(config['DE_NOVO'])

    if config['TRAIN'] == True: # Train

        db = dataload(config['DB'])
        print('>> Dataload done ...\n')

        de_novo = denovo_parsing(de_novo)
        db = db_parsing(db)
        clu = new_clu(cluster_info, mgf_info, config['PPM'])
        df = merge(de_novo, clu, db)

    else: # Test

        print('>> Dataload done ...\n')

        de_novo = denovo_parsing(de_novo)
        clu = new_clu(cluster_info, mgf_info, config['PPM'])
        df = merge(de_novo, clu)

    print('>> Generating candidates ...')

    # df_10 = create_dataset(df, 10, config['TRAIN'])
    df_2 = create_dataset(df, 2, config['TRAIN'])

    print('>> Done \n')

    # feature extraction
    print('>> Feature extraction ...')

    df_feature = feature(df_2)
    df_feature = ifi_feature(df_feature, remove_path)
    df_feature['RT'] = df_feature['RT'] / (config['ELUTION_TIME']/60)
    df_feature = rt_feature_training(df_feature)

    df_feature.to_csv(config['RESULT_NAME'], index=False)

    print()

    cometX_mgf_path = config['MGF_XCORR']
    if os.path.exists(cometX_mgf_path):
        shutil.rmtree(cometX_mgf_path)
    os.mkdir(cometX_mgf_path)

    if config['TRAIN'] == True: # Train
        train_df = df_feature[df_feature['GT'].notnull()]
        CometX_mgf(train_df, config['MGF_PATH'], cometX_mgf_path)
    else:
        CometX_mgf(df_feature, config['MGF_PATH'], cometX_mgf_path)

    print()
    print('>> Done ...\n')

    # Calculatie XCorr
    print('We need to calculate the XCorr value')