import pandas as pd

def File_index(mgf_list):
    
    dic = {}
    
    for idx, i in enumerate(mgf_list):
        dic[idx] = i
        
    return dic

def Cluster_csv(path_dir, cluster_list, index_dic):
    
    cnt = 0
    cluster = []
    append = cluster.append

    for f in cluster_list:

        clu = open(path_dir+'\\'+f).readlines()

        for i in clu:
            if cnt == 0:
                clu_cnt = int(i.split('\t')[1])
                clu_num = i.split('\t')[0].split('.')[1]
                cnt += 1
            elif i == '\n':
                cnt = 0
            else:
                idx = int(i.split('\t')[1])
                fn = index_dic[idx]
                scan = int(i.split('\t')[2])

                temp = fn, scan, clu_cnt, clu_num

                append(temp)
                
    return pd.DataFrame(cluster, columns=['Source File','Scan number','count','cluster']).sort_values(['Source File','Scan number']).reset_index(drop=True)

