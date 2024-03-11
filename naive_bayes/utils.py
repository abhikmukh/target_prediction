from itertools import islice
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Any, Iterator, Iterable, List, Tuple
from functools import lru_cache
import os
import csv

from rdkit import Chem
from rdkit.Chem import AllChem


def batch_iterable(iterable: Iterable[Any], batch_size: int) -> Iterator[List[Any]]:
    iterator = iter(iterable)
    while batch := list(islice(iterator, batch_size)):
        yield batch


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


@lru_cache(maxsize=10000)
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


def get_fp_chunk_list(chunk: List[Tuple[str, str, str]]) -> List[Tuple[str, str, str, list]]:
    result: List[Tuple[str, str, str, list]] = []

    for row in chunk:
        try:
            fp = get_fp_rdkit(row[1])
            result.append((row[0], row[1], row[2], list(fp)))
        except:
            pass
    return result
