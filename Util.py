# Importing Modules
import numpy as np
import os

def match(dataset):
    print('>> Start ...')
    print()

    match = 0
    shift = True

    total = len(dataset)

    for idx, (a, b, i, j) in enumerate(dataset[['Source File', 'Scan number', 'Sequence', 'GT']].values):

        if idx == 0:
            key = a, b

        if key == (a, b):
            if shift == True:
                if i == j:
                    match += 1
                    shift = False
        else:
            key = a, b
            shift = True
            if i == j:
                match += 1
                shift = False

    print('>> Matching PSM :', match)

def Global_Alignment(sequence_1, sequence_2):
    # Needleman-Wunsch Alignment

    # Creat Matrices
    main_matrix = np.zeros((len(sequence_1) + 1, len(sequence_2) + 1))
    match_checker_matrix = np.zeros((len(sequence_1), len(sequence_2)))

    # Provideing the scores for match, mismatch and gap
    match_reward = 1
    mismatch_penalty = -1
    gap_penalty = -2

    # Fill the match ckecker matrix accorording to match or mismatch
    for i in range(len(sequence_1)):
        for j in range(len(sequence_2)):
            if sequence_1[i] == sequence_2[j]:
                match_checker_matrix[i][j] = match_reward
            else:
                match_checker_matrix[i][j] = mismatch_penalty
    # print(match_checker_matrix)

    # Filling up the matrix using Needleman_Wunsch algorithm
    # Step 1 : Initialization

    for i in range(len(sequence_1) + 1):
        main_matrix[i][0] = i * gap_penalty
    for j in range(len(sequence_2) + 1):
        main_matrix[0][j] = j * gap_penalty

    # Step 2 : Matrix Filling

    for i in range(1, len(sequence_1) + 1):
        for j in range(1, len(sequence_2) + 1):
            main_matrix[i][j] = max(main_matrix[i - 1][j - 1] + match_checker_matrix[i - 1][j - 1],
                                    main_matrix[i - 1][j] + gap_penalty,
                                    main_matrix[i][j - 1] + gap_penalty)
    #     print(main_matrix)

    # Step 3 :  Traceback

    aligned_1 = ''
    aligned_2 = ''

    ti = len(sequence_1)
    tj = len(sequence_2)

    # print(ti,tj)

    while (ti > 0 or tj > 0):

        if (ti > 0 and tj > 0 and main_matrix[ti][tj] == main_matrix[ti - 1][tj - 1] + match_checker_matrix[ti - 1][tj - 1]):
            aligned_1 = sequence_1[ti - 1] + aligned_1
            aligned_2 = sequence_2[tj - 1] + aligned_2

            ti = ti - 1
            tj = tj - 1

        elif (ti > 0 and main_matrix[ti][tj] == main_matrix[ti - 1][tj] + gap_penalty):
            aligned_1 = sequence_1[ti - 1] + aligned_1
            aligned_2 = '-' + aligned_2

            ti = ti - 1

        else:
            aligned_1 = '-' + aligned_1
            aligned_2 = sequence_2[tj - 1] + aligned_2

            tj = tj - 1

    #     print(ti,tj)

    # test

    #     print(aligned_1)
    #     print(aligned_2)
    return aligned_1, aligned_2


def Match_score(seq_1, seq_2):
    return list(np.array(list(seq_1)) == np.array(list(seq_2))).count(True)


def top1(dataset):
    dataset = dataset.drop_duplicates(subset=['Source File', 'Scan number'], keep='first').reset_index(drop=True)

    return dataset

def amino_acid_recall(dataset, gt='GT', pep='Peptide'):

    cnt = 0
    for i in dataset[[gt]].values:
        cnt += len(i[0])

    top_1 = top1(dataset)
    top_1 = top_1[top_1[pep].notnull()]

    if 'Sequence' in top_1.columns:
        col = 'Sequence'
    else:
        col = pep

    aa_recall = 0
    for idx, (i, j) in enumerate(top_1[[col, gt]].values):

        if idx % 100000 == 0:
            print(idx)

        if type(i) is not str or type(j) is not str:
            pass
        else:
            a, b = Global_Alignment(i.replace('m', 'M'), j)
            aa_recall += Match_score(a, b)

    print(round(aa_recall/cnt, 4))

def Make_txt(mgf_path, save_path, file_name):

    mgf_list = os.listdir(mgf_path)

    f = open(save_path + "\\" + file_name, 'w')

    for i in mgf_list:
        f.write(mgf_path + '\\' + i + '\n')
    f.close()