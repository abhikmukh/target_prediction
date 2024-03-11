#!/usr/bin/env python
# coding: utf-8

# In[17]:


import psycopg2
import psycopg2.extras
import math
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import multiprocessing
import csv
import time




def read_csv(filename):
    '''
    ead csv file
    :param filename:
    :return:
    '''
    with open(filename, 'r') as csv_file:
        next(csv_file)
        csv_reader = csv.reader(csv_file)
        for rows in csv_reader:
            yield rows




def get_fp_rdkit(smiles: str) -> set:
    """Function to calculate MorganFingerprint from smiles.
    It returns index of all '1' bits of not-folded fingerprint.
    Args:
        smiles (str): smiles string
    Returns:
        set: return list of index of '1' bits.
    """

    mol = Chem.MolFromSmiles(smiles)

    if not mol:
        return

    fp = AllChem.GetMorganFingerprint(mol, 2)
    if not fp:
        print(smiles)
        return

    return set(fp.GetNonzeroElements().keys())


def calculate_correction_np(target_id, fingerprint_list, freq_all_list, fingerprint_freq_all_list):
    '''

    Args:
      target_id: #

    Returns:

    '''

    conn = psycopg2.connect("host=localhost dbname=nb_target user=postgres password=abhik1234")

    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    cur.execute("""select 
    fp.feature as present
    from 
    normalised_probabilities np,
    fingerprints fp,
    targets pt
    where 
    np.target_id = pt.target_id
    and
    pt.sc_target_id = '{target_id}'
    and
    np.fingerprint_id = fp.fingerprint_id
    and
    fp.feature in {ids}
    """.format(target_id=target_id, ids=tuple(fingerprint_list)))

    present_fp = cur.fetchall()
    present_fp_list = []
    for x_fp in present_fp:
        present_fp_list.append(x_fp[0])

    cur.execute("""select
    pt.laplaciank
    from targets as pt
    where
    pt.sc_target_id = '{target_id}'
    """.format(target_id=target_id))

    laplacian = cur.fetchall()
    np_corr_list = []
    for fp, freq_all in zip(fingerprint_freq_all_list, freq_all_list):

        if fp not in present_fp_list:
            denominator = ((freq_all + laplacian[0][0])*(1/laplacian[0][0]))

            correction = math.log(1/denominator)

            np_corr_list.append(correction)

    return sum(np_corr_list)


def predict_target(smiles):
    conn = psycopg2.connect("host=localhost dbname=nb_target user=postgres password=abhik1234")

    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    if get_fp_rdkit(smiles) is not None:

        fingerprint_list = list(get_fp_rdkit(smiles))

        # 1. get sum(NP) for known fp
        cur.execute("""select * from (select pt.sc_target_id, 
        sum(np.normalised_probability) as np_sum
        from 
        normalised_probabilities np,
        fingerprints fp,
        targets pt
        where 
        np.fingerprint_id = fp.fingerprint_id
        and
        np.target_id = pt.target_id
        and
        fp.feature in {ids}
        group by 
        pt.sc_target_id
        order by
        np_sum desc) as a
        fetch first 100 rows only
        """.format(ids=tuple(fingerprint_list)))

        target_id_score = cur.fetchall()
        np_score_dict_list = []
        for row in target_id_score:
            np_score_dict_list.append(dict(row))


        target_id_list = []
        np_score_list = []
        for my_dict in np_score_dict_list:
            target_id_list.append(my_dict['sc_target_id'])
            np_score_list.append(my_dict['np_sum'])

        # 2. Get array of known FP


        cur.execute("""select 
            fp.feature,
            fp.frequency_all
            from 
            fingerprints fp
            where 
            fp.feature in {ids}
            """.format(ids=tuple(fingerprint_list)))

        frequency_all_array = cur.fetchall()

        freq_all_list = []
        fingerprint_freq_all_list = []
        for x in frequency_all_array:
            freq_all_list.append(x[1])
            fingerprint_freq_all_list.append((x[0]))




        target_id_score = []

        for target_id, np_score in zip(target_id_list, np_score_list):
            np_corr_per_target = (calculate_correction_np(target_id, fingerprint_list, freq_all_list=freq_all_list, fingerprint_freq_all_list=fingerprint_freq_all_list))

            tup = (target_id, float(np_score), np_corr_per_target, (float(np_score) + np_corr_per_target))
            target_id_score.append(tup)

        df = pd.DataFrame(target_id_score, columns=["target_id", "np_score", "correction", "np_score_corrected_sum"])
        sorted_df = (df.sort_values(by=['np_score_corrected_sum'], ascending=False))
        result_df = (sorted_df.drop(columns=["np_score", "correction", "np_score_corrected_sum"]))
        target_list = (result_df["target_id"].to_list()[:10])
        return smiles, target_list
    else:
        print("failed" + "_" + smiles)


def main():
    start = time.time()
    p = multiprocessing.Pool(8)
    data = read_csv("test_data.csv")
    #result_dict = {}
    smiles_list = []
    target_list = []
    for x in data:
        print(x)
        smiles = (x[1])
        smiles_list.append(smiles)
        target = x[0]
        target_list.append(target)
    result_list = p.map(predict_target, smiles_list)
    print(len(result_list))
    df = pd.DataFrame(result_list)
    df.to_csv("prediction_list_poc-2.csv", index=False)

    end = time.time()
    print(end - start)

    #
    # for result in result_list:
    #
    #     for target, smiles in zip(target_list, smiles_list):
    #         print(smiles)
    #         print(result[0])
    #
    #         if result:
    #             if smiles == result[0]:
    #                 if target in result[1]:
    #                     rank = result[1].index(target)
    #                 else:
    #                     rank = 5
    #                 result_dict[target] = rank
    #
    # print(len(result_dict))
    # print(result_dict)
    #
    # with open('target_prediction.csv', 'w', newline='') as f:
    #     # Create a CSV writer object
    #     writer = csv.writer(f)
    #     # Write one key-value tuple per row
    #     for row in result_dict.items():
    #         writer.writerow(row)


if __name__ == "__main__":
    main()
