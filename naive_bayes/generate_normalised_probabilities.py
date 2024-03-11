import pandas as pd
import os
import math
from collections import defaultdict, Counter

from .utils import read_csv


def write_laplacian_df(num_of_feature_in_a_class, total_num_mols, output_file):

    p_baseline_dict_for_all_classes = {k:v / total_num_mols for k, v in num_of_feature_in_a_class.items()}

    laplacian_df = pd.DataFrame(p_baseline_dict_for_all_classes.items(), columns=['sc_target_id', 'baseline'])

    laplacian_df['laplaciank'] = laplacian_df['baseline'].apply(lambda x: 1 / x)
    laplacian_df = laplacian_df.drop(columns=['baseline'])
    print(laplacian_df.head())
    laplacian_df.to_csv(output_file, index=False)


def create_total_feature_counts_dict(target_fingerprint_csv_file):
    rows = read_csv(target_fingerprint_csv_file)
    feature_list = []
    for row in rows:
        feature_list.append(int(row[1]))
    total_feature_counts_dict = Counter(feature_list)
    return total_feature_counts_dict


def write_map_fingerprint(total_feature_counts_dict, output_file):

    total_feature_count_df = pd.DataFrame(total_feature_counts_dict.items(), columns=["feature", "frequency_all"])
    total_feature_count_df.to_csv(output_file, index=False, chunksize=100000)


def create_result_dict(target_fingerprint_csv_file):
    result_dict = defaultdict(lambda: defaultdict(int))
    rows = read_csv(target_fingerprint_csv_file)
    for row in rows:
        result_dict[row[0]][row[1]] += 1
    return result_dict


def write_target_norm_p(result_dict, num_of_feature_in_a_class, len_concat_list, total_feature_counts_dict,
                        output_file):
    list_of_tuple = []
    logp_final_new = defaultdict(dict)

    for k, v in result_dict.items():

        for x, y in v.items():


            if num_of_feature_in_a_class.get(k):
                p_base = (num_of_feature_in_a_class.get(k)) / len_concat_list
                numerator = (y + 1)

                denominator = (total_feature_counts_dict.get(int(x)) + (1 / p_base))

                p_corr = numerator / denominator
                norm_p = math.log(p_corr / p_base)
                logp_final_new[k][x] = norm_p
                tup = (k, x, norm_p)
                list_of_tuple.append(tup)
    target_norm_df = pd.DataFrame(list_of_tuple, columns=['target', 'fingerprint', 'norm_prob'])
    target_norm_df.to_csv(output_file, index=False)


def run(data_dir: str = "data") -> None:
    col_names = ["monomer_id", "smiles", "target", "sets"]
    df = pd.read_csv(
        (os.path.join(data_dir, "train_fp.csv")),
        names=col_names)

    df.to_csv((os.path.join(data_dir, "map_compounds_targets.csv")), columns=["smiles", "target"], index=False)
    df.to_csv((os.path.join(data_dir, "map_compounds.csv")), columns=["monomer_id", "smiles"], index=False)

    total_feature_counts_dict = create_total_feature_counts_dict(
        target_fingerprint_csv_file=(os.path.join(data_dir, "target_fp.csv")))

    write_map_fingerprint(total_feature_counts_dict, output_file=(os.path.join(data_dir, "map_fingerprints.csv")))

    target_fp_df_col_names = ["target", "fingerprint"]

    target_fp_df = pd.read_csv(
        (os.path.join(data_dir, "train_fp.csv")), names=target_fp_df_col_names)
    target_fp_shape = target_fp_df.shape
    len_concat_list = target_fp_shape[0]

    num_of_feature_in_a_class = dict(target_fp_df.groupby('target').fingerprint.count())
    write_laplacian_df(num_of_feature_in_a_class=num_of_feature_in_a_class, total_num_mols=len_concat_list,
                       output_file=(os.path.join(data_dir, "map_target_laplacian.csv")))

    result_dict = create_result_dict(
        target_fingerprint_csv_file=(os.path.join(data_dir, "target_fp.csv")))

    write_target_norm_p(result_dict, num_of_feature_in_a_class, len_concat_list, total_feature_counts_dict,
                        output_file=(os.path.join(data_dir, "map_target_norm_p.csv")))


if __name__ == "__main__":
    run()
