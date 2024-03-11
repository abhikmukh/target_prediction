import os

import pandas as pd
from sklearn.model_selection import train_test_split


class TrainTestPreparation:
    """
    This class is used to create training and test set from the given csv file  and save them in the data directory

    """
    def __init__(self, data_dir: str) -> None:
        self.data_dir = data_dir

    def create_training_and_test_set(self, csv_file: str) -> None:
        csv_file_path = os.path.join(self.data_dir, csv_file)
        df = pd.read_csv(csv_file_path)
        train, test = train_test_split(df, test_size=0.1, random_state=0, stratify=df[['target']])
        print(f"Shape of the training set is {train.shape}")
        print(f"Shape of the test set is {test.shape}")
        train.to_csv(os.path.join(self.data_dir, "train.csv"), index=False)
        test.to_csv(os.path.join(self.data_dir, "test.csv"), index=False)
        test_modified = test.drop(columns=["target"])
        test_modified.drop_duplicates(inplace=True)
        print(f"Number of unique compounds in the test set is {test_modified.shape[0]}")
        test_modified.to_csv(os.path.join(self.data_dir, "test_set_unique_compounds.csv"), index=False)
        print("Training and test set created successfully")


