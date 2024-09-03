import pandas as pd
import numpy as np
import re

from sklearn.cluster import DBSCAN

def dataload(file_name):

    f_e = file_name.split('.')[-1]
    
    if f_e == 'csv':
        chunk_data = pd.read_csv(file_name, chunksize=1000000)
        data = pd.concat([chunk for chunk in chunk_data])

    if 'Scan number' in data.columns:
        data = data.astype({'Scan number': 'int'})

    return data

def top1(dataset):
    
    dataset = dataset.drop_duplicates(subset=['Source File', 'Scan number'], keep='first').reset_index(drop=True)
    
    return dataset

def strip_seq(seq):
    
    # pep = ''.join([_ for _ in seq if ord(_) in range(65, 91)]) # 65 ~ 91 대문자
    pep = re.compile('[^0-9()+-.]')
    
    res = ''.join(pep.findall(seq)).replace('I', 'L')
    res = res.replace('m', 'M')
    
    return res

def strip_sequence(de_novo, col, col_):
    
    de_novo[col_] = de_novo[col].apply(strip_seq)
    
    return de_novo


def denovo_parsing(de_novo):

    print('>> De novo parsing strat ...')

    de_novo = de_novo[de_novo['Peptide'].notnull()] # Double Check

    ''' 
        1. using strip sequence
        2. AA I equal to L
                
    '''

    de_novo = strip_sequence(de_novo, 'Peptide', 'Sequence')
    
    print('>> De novo parsing done \n')
    
    return de_novo.sort_values(by=['Source File', 'Scan number']).reset_index(drop=True) # Double Check

def same_mass(seq):
    
    return seq.replace('I', 'L')

def i_to_l(db, col, col_):
    
    db[col_] = db[col].apply(same_mass)
        
    return db

def db_parsing(db):
    
    print('>> DB parsing strat ...')
    
    ''' 
        1. AA I equal to L
                
    '''
        
    db = i_to_l(db, 'GT', 'GT')
    db.sort_values(by=['Source File', 'Scan number']).reset_index(drop=True)
    print('>> DB parsing done \n')
    
    return db.sort_values(by=['Source File', 'Scan number']).reset_index(drop=True)

def remove_no_clu_info(dataset):
        
    return dataset.dropna(subset=['cluster']).reset_index(drop=True)

def clu_plus_charge(dataset):
        
    print('>> Cluster number + Charge ...')
    
    temp = [] 
    append = temp.append
        
    for idx, (i, j) in enumerate(dataset[['cluster', 'z']].values):
        
        append(str(i)+'_'+str(int(j)))
    
    return temp

def clu_num_to_int(clu):
    
    print('>> Change to integer value ...')
    
    ind = 0
    dic = {}
    
    temp2 = set(clu)
    
    for i in temp2:
        if i not in dic:
            dic[i] = ind
            ind += 1
    
    # integer cluster number
    
    temp = [] 
    append = temp.append
    
    for i in clu:
        append(dic[i])
    
    return temp

def change_clu_num(dataset):
   
    print('>> Cluster index pre-processing ...')

    '''
        MScluster result does not consider the charge value
        Seperate the cluster by adding charge
    '''
    
    temp = clu_plus_charge(dataset)

    
    ''' 
        Change the modified cluster number (cluster number + charge) to integer value
    '''
    
    dataset['cluster'] = clu_num_to_int(temp)
    
    return dataset

def merge(de_novo, cluster, db = None):
    
    if db is not None:
        print('>> De novo + Database searching result merging ...')
        print('>> Making Training, Test set ...')
        df = pd.merge(de_novo, db, how='outer')
        
        print('>> Clustering result merging ... \n')
        df = pd.merge(df, cluster, how='left')
        
        top_1 = top1(df)
        print('>> Reliable PSM (GT) :', top_1['GT'].notnull().sum(), 'Scans \n')
        
        # Clustering 정보 없는 것 고려대상에서 제외, 제거
        df = remove_no_clu_info(df)

    else:
        print('>> De novo + Clustering result merging ... \n')
        df = pd.merge(de_novo, cluster, how='left')
        df = remove_no_clu_info(df)

    print('>> Merging done \n')

    return df.reset_index(drop=True)

# Original Cluster Size : count
# New Cluster Size : new_count

def refinement(dataset, col, eps, min_samples):
        
    feature = dataset[[col, 'm/z']] # 'RT'

    model = DBSCAN(eps=eps, min_samples=min_samples, n_jobs = -1)
    predict = model.fit_predict(feature)
    
    dataset['new_clu'] = predict
    
    return dataset

def adding_top10(top1, top10, col):
    
    top_1 = top1[['Source File', 'Scan number', col]]

    merged_dataset = pd.merge(top10, top_1)
        
    return merged_dataset

def clu_size(top1, col):
    
    dic = {}
    clu_list = top1[col]
    
    for i in clu_list:
        if i not in dic:
            dic[i] = 0
        dic[i] += 1
        
    temp = [] 
    append = temp.append
    
    for i in clu_list:
        append(dic[i])
    
    top1['new_count'] = temp
    
    return top1

def cluster_refinement(dataset, col, eps, min_samples):
    
    print('>> Refinement ...')
    
    top_1 = top1(dataset)
    top_1 = refinement(top_1, col, eps, min_samples)
    
    dataset = adding_top10(top_1, dataset, 'new_clu')
    
    top_1 = clu_size(top_1, 'new_clu')
    dataset = adding_top10(top_1, dataset, 'new_count')
    
    print('>> Refinement done \n')
    
    return dataset

def new_clu(cluster, info, ppm):
    
    temp = info['m/z']*1e-6*ppm
    value = np.median(list(temp))
    
    df = pd.merge(info, cluster, how='left')
    
    df = change_clu_num(df)
    df = cluster_refinement(df, 'cluster', value, 1)
    
    return df

def top_n(dataset, n):
    
    dic = {}
    
    dataset['Sequence'] = dataset['Sequence'].fillna('nan')
    dataset['Score'] = dataset['Score'].fillna(0)
    
    for i, j, k in dataset[['new_clu', 'Sequence', 'Score']].values:
            
        if i not in dic:
            dic_temp = {}
            dic[i] = dic_temp

        if j not in dic[i]:
            dic[i][j] = 0
        dic[i][j] += k*0.01
    
    topn = {}

    for i in dic:
        temp = sorted(dic[i].items(), reverse=True, key=lambda x: x[1])
        temp = temp[0:n]

        for idx, j in enumerate(temp):
            if i not in topn:
                topn[i] = []
            topn[i].append((j, idx))
    
    return topn

def candidate_info(dataset, candi):
    
    dic_pep, dic_seq, dic_rank, dic_score = {}, {}, {}, {}

    for i, j, k in dataset[['Peptide', 'Sequence', 'new_clu']].values:
        for t in candi[k]:
            if j in t[0][0]:
                if k not in dic_pep:
                    dic_seq[k] = []
                    dic_pep[k] = []
                    dic_rank[k] = []
                    dic_score[k] = []

                if i not in dic_pep[k]:
                    dic_seq[k].append(j)
                    dic_pep[k].append(i)
                    dic_rank[k].append(t[1])
                    dic_score[k].append(t[0][1])
    
    return dic_pep, dic_seq, dic_rank, dic_score

def candidates(dataset, dic_p, dic_seq, dic_rank, dic_score):

    temp = []
    append = temp.append

    for i, j, k in dataset[['Source File', 'Scan number', 'new_clu']].values:
        for l, t, z, y in zip(dic_p[k], dic_seq[k], dic_rank[k], dic_score[k]):
            data = i, j, l, t, y, z
            append(data)
    
    return pd.DataFrame(temp, columns=['Source File', 'Scan number', 'Peptide', 'Sequence', 'Score', 'Rank'])

def labeling(dataset, t=0, f=1): # seq=gt : 0, otherwise : 1

    temp = []
    append = temp.append

    for i, j in dataset[['Sequence', 'GT']].values:
        if i == j:
            append(t)
        else:
            append(f)
    
    dataset['Label'] = temp
    
    return dataset

def remove_psm(dataset, train):

    temp1 = dataset['Peptide'].notnull()

    if train == True:
        temp2 = dataset['GT'].notnull()
        return dataset[temp1|temp2]

    return dataset[temp1]


def create_dataset(dataset, n, train):
    
    n_candidates = top_n(dataset, n)
    pep, seq, rank, score = candidate_info(dataset, n_candidates)
    
    df_1 = top1(dataset)
    dataset = candidates(df_1, pep, seq, rank, score)

    if train == True:
        df_1_ = df_1[['Source File', 'Scan number', 'GT', 'z', 'm/z', 'RT', 'new_count']]
        dataset = pd.merge(dataset, df_1_, how='outer')
        dataset = labeling(dataset)
    else:
        df_1_ = df_1[['Source File', 'Scan number', 'z', 'm/z', 'RT', 'new_count']]
        dataset = pd.merge(dataset, df_1_, how='outer')

    dataset = remove_psm(dataset, train)

    dataset = dataset.sort_values(by=['Source File', 'Scan number', 'Rank']).reset_index(drop=True) # Double Check
    
    return dataset

def delta_score(dataset):
    
    temp = []
    append = temp.append

    for idx, (i, j, z) in enumerate(dataset[['Source File', 'Scan number', 'Score']].values):

        if idx == 0:
            key = (i, j)
            append(0)
            cri = z
            continue

        if key == (i, j):
            append(cri-z)
        else:
            key = (i, j)
            append(0)
            cri = z
                
    return temp

def log_scale(dataset):
    
    dataset['Score'] = np.log(dataset['Score']+1)
    dataset['new_count'] = np.log(dataset['new_count']+1)
    
    return dataset

def feature(dataset):
    
    dataset = log_scale(dataset)
    dataset['delta'] = delta_score(dataset)

    return dataset.reset_index(drop=True)