from itertools import islice
from multiprocessing import Pool, cpu_count
from typing import Any, Iterator, Iterable, List, Tuple
import os
import csv


from .utils import batch_iterable, read_csv, get_fp_chunk_list


class FingerprintGenerator:
    def __init__(self, data_dir):
        self.data_dir = data_dir
        self.cpu_count = cpu_count()

    def _get_batched_smiles(self, csv_file):
        chunk_size = 1024
        data_file = os.path.join(self.data_dir, csv_file)
        smiles_data = read_csv(data_file)
        return batch_iterable(smiles_data, chunk_size)

    def generate_fingerprint_set(self, csv_file, output_file):

        output_file = os.path.join(self.data_dir, output_file)
        batched_smiles = self._get_batched_smiles(csv_file)
        with open(output_file, "w") as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            with Pool(processes=self.cpu_count) as pool:
                for results in pool.imap_unordered(get_fp_chunk_list, batched_smiles):
                    for result in results:
                        writer.writerow(result)

    def generate_fingerprint_set_new(self, csv_file, output_file):

        output_file = os.path.join(self.data_dir, output_file)
        batched_smiles = self._get_batched_smiles(csv_file)
        with open(output_file, "w") as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            with Pool(processes=self.cpu_count) as pool:
                for results in pool.imap_unordered(get_fp_chunk_list, batched_smiles):

                    for result in results:
                        # result = ('192083', 'COC(=O)Cn1nc(-c2cccc(OCCC(C)(C)C)c2)c2ccccc12',
                        # 'Cytochrome b', [2154935424, ...])
                        if result[3]:
                            for fingerprint in result[3]:
                                # csv_file.write(f"'{result[2][0]}', {fingerprint}\n")
                                if "," in result[2]:
                                    target = result[2].split(",")[0]
                                    csv_file.write(f"{target}, {fingerprint}\n")
                                else:
                                    csv_file.write(f"{result[2]}, {fingerprint}\n")






