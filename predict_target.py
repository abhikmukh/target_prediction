import psycopg2
import psycopg2.extras
import math
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import multiprocessing
import csv
import time

from naive_bayes.utils import  read_csv, get_fp_rdkit
from naive_bayes.database_utils import DatabaseConnector

base_path = os.getcwd()
database_connector = DatabaseConnector(base_path=base_path, cred_file="db_credential.yaml")


def _calculate_correction_np(target_id, fingerprint_list, freq_all_list, fingerprint_freq_all_list):

    query_to_get_fps = f"select fp.feature as present from normalised_probabilities np, fingerprints fp, targets pt " \
            f"where np.target_id = pt.target_id and pt.sc_target_id = '{target_id}' " \
            f"and np.fingerprint_id = fp.fingerprint_id and fp.feature in {tuple(fingerprint_list)}"

    present_fp = database_connector.run_query(query_to_get_fps)
    present_fp_list = []
    for fp in present_fp:
        present_fp_list.append(fp[0])

    query_to_get_laplacian = f"select pt.laplaciank from targets as pt where pt.sc_target_id = '{target_id}'"
    laplacian = database_connector.run_query(query_to_get_laplacian)

    np_corr_list = []
    for fp, freq_all in zip(fingerprint_freq_all_list, freq_all_list):

        if fp not in present_fp_list:
            denominator = ((freq_all + laplacian[0][0])*(1/laplacian[0][0]))

            correction = math.log(1/denominator)

            np_corr_list.append(correction)

    return sum(np_corr_list)


def predict_target(input_smiles):

    if get_fp_rdkit(input_smiles) is not None:

        fingerprint_list = list(get_fp_rdkit(input_smiles))

        # 1. get sum(NP) for known fp
        query_target_id_score = f"select * from (select pt.sc_target_id, sum(np.normalised_probability) as np_sum " \
                                f"from normalised_probabilities np, fingerprints fp, targets pt " \
                                f"where np.fingerprint_id = fp.fingerprint_id and np.target_id = pt.target_id " \
                                f"and fp.feature in {tuple(fingerprint_list)} group by pt.sc_target_id " \
                                f"order by np_sum desc) as a " \
                                f"fetch first 100 rows only"

        target_id_score = database_connector.run_query(query_target_id_score)

        results_as_tuples = [(str(row[0]), float(row[1])) for row in target_id_score]

        np_score_dict_list = [dict(zip(('sc_target_id', 'np_sum'), row)) for row in results_as_tuples]

        target_id_list = []
        np_score_list = []
        for my_dict in np_score_dict_list:
            target_id_list.append(my_dict["sc_target_id"])
            np_score_list.append(my_dict["np_sum"])

        # 2. Get array of known FP

        query_frequency_all_array = f"select fp.feature, fp.frequency_all from fingerprints fp " \
                                    f"where fp.feature in {tuple(fingerprint_list)}"
        frequency_all_array = database_connector.run_query(query_frequency_all_array)

        freq_all_list = []
        fingerprint_freq_all_list = []
        for x in frequency_all_array:
            freq_all_list.append(x[1])
            fingerprint_freq_all_list.append((x[0]))

        target_id_score = []

        for target_id, np_score in zip(target_id_list, np_score_list):
            np_corr_per_target = (_calculate_correction_np(target_id=target_id, fingerprint_list=fingerprint_list,
                                                           freq_all_list=freq_all_list,
                                                           fingerprint_freq_all_list=fingerprint_freq_all_list))

            tup = (target_id, float(np_score), np_corr_per_target, (float(np_score) + np_corr_per_target))
            target_id_score.append(tup)

        df = pd.DataFrame(target_id_score, columns=["target_id", "np_score", "correction", "np_score_corrected_sum"])
        sorted_df = (df.sort_values(by=['np_score_corrected_sum'], ascending=False))
        result_df = (sorted_df.drop(columns=["np_score", "correction", "np_score_corrected_sum"]))
        target_list = (result_df["target_id"].to_list()[:5])
        return str(target_list)
    else:
        print("failed" + "_" + input_smiles)




