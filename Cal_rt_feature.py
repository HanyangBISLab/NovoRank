import os
import warnings
import numpy as np
import pandas as pd

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

warnings.filterwarnings(action='ignore')

from tqdm import tqdm
from deeplc import DeepLC

def top2(dataset):
    
    dataset = dataset.sort_values(by=['Source File', 'Scan number', 'Rank'])
    
    dataset1 = dataset.drop_duplicates(subset=['Source File', 'Scan number'], keep='first')
    dataset2 = dataset.drop_duplicates(subset=['Source File', 'Scan number'], keep='last')
    
    new = pd.concat([dataset1,dataset2]).sort_values(by=['Source File', 'Scan number', 'Rank']).reset_index(drop=True)
    
    return new

def cal_cal(dataset):
    
    dataset = dataset.drop_duplicates(subset=['Source File', 'Scan number'], keep='first')
    dataset = dataset.drop_duplicates(subset=['seq', 'modifications'], keep='first')
    dataset = dataset.sort_values(by=['Score'], ascending=False).reset_index(drop=True)
    dataset = dataset.iloc[:1000]
    dataset = dataset[['seq', 'modifications', 'tr']]
    
    return dataset

def pred_rt_data(dataset):
    
    mod = []
    append = mod.append
    
    dataset = dataset[['Source File', 'Scan number', 'Sequence', 'Peptide', 'RT', 'Score']]
    
    for i in dataset[['Peptide']].values:
        mo = ""
        for idx, j in enumerate(i[0]):
            if j == "C":
                mo=mo+str(idx+1)+'|'+"Carbamidomethyl"+'|'
            elif j == 'm':
                mo=mo+str(idx+1)+'|'+"Oxidation"+'|'

        append(mo[:len(mo)-1])
        
    dataset['modifications'] = mod
    dataset.columns =['Source File', 'Scan number', 'seq', 'Peptide', 'tr', 'Score', 'modifications']
    
    cal = dataset[['Source File', 'Scan number', 'seq', 'modifications', 'tr', 'Score']]
    cal['modifications'] = ""
    
    pep = dataset[['Source File', 'seq', 'modifications', 'tr']]
    
    cal['modifications'] = cal['modifications'].fillna("")
    pep['modifications'] = pep['modifications'].fillna("")
    
    return pep, cal

def pred_rt_cal(pep, cal):

    pred_rt = []

    for i in tqdm(sorted(set(pep['Source File']))):

        pep_ = pep[pep['Source File'] == i]
        cal_ = cal[cal['Source File'] == i]
        
        pep_ = pep_[['seq', 'modifications', 'tr']]
        cal_ = cal_cal(cal_)
        
        dlc = DeepLC()
        
        dlc.calibrate_preds(seq_df=cal_)
        preds = dlc.make_preds(seq_df=pep_)

        pred_rt += preds
    
    return pred_rt

def rt_feature_training(dataset):
    
    dataset_ = dataset[dataset['Peptide'].isnull()]
    dataset = dataset[dataset['Peptide'].notnull()]
    dataset = top2(dataset)
    
    pep, cal = pred_rt_data(dataset)
    res = pred_rt_cal(pep, cal)
    
    dataset['pred_RT'] = res
    dataset['diff_RT'] = abs(dataset['RT']-dataset['pred_RT'])
    dataset['diff_RT'] = np.log(dataset['diff_RT']+1)
    
    dataset_['pred_RT'] = np.nan
    dataset_['diff_RT'] = np.nan
    
    dataset = pd.concat([dataset, dataset_]).sort_values(by=['Source File', 'Scan number', 'Rank']).reset_index(drop=True)

    return dataset