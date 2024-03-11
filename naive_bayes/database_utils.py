import yaml
import os
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy import text
from typing import Sequence

class DatabaseConnector:
    def __init__(self, base_path: str, cred_file: str) -> None:
        self.base_path = base_path
        self.cred_file = cred_file

    def read_db_creds(self):
        """
        Read database credentials from yaml file
        """
        full_file_path = os.path.join(self.base_path, self.cred_file)
        with open(full_file_path, 'r') as file:
            data = yaml.safe_load(file)
            return data

    def init_db_engine(self):
        """
        Initialize database engine
        """
        database_cred = self.read_db_creds()
        database_uri = f"postgresql+psycopg2://{database_cred['USER']}" \
                       f":{database_cred['PASSWORD']}@{database_cred['HOST']}:" \
                       f"{database_cred['PORT']}/{database_cred['DATABASE']}"
        database_engine = create_engine(database_uri)
        return database_engine

    def upload_to_db(self, csv_file, table_name):
        """
        Upload data from csv file to database table
        """
        engine = self.init_db_engine()
        df = pd.read_csv(csv_file)

        df.to_sql(name=table_name, con=engine, if_exists='replace', index=False)
        print(f"Data uploaded to table {table_name}")

    def alter_table_data_type(self, table_name: str, column_name: str, data_type: str) -> None:
        """
        Alter table column data type
        """
        engine = self.init_db_engine()
        with engine.connect() as con:
            con.execute(text(f"ALTER TABLE {table_name} ALTER COLUMN {column_name} TYPE {data_type};"))
            con.commit()
        print(f"Column {column_name} of type {data_type} updated to table {table_name}")

    def add_column_with_primary_key(self, table_name: str, column_name: str, data_type: str) -> None:
        """
        Add column to table with primary key
        """
        engine = self.init_db_engine()
        with engine.connect() as con:
            con.execute(text(f"ALTER TABLE {table_name} ADD COLUMN {column_name} {data_type} PRIMARY KEY;"))
            con.commit()
        print(f"Column {column_name} of type {data_type} added to table {table_name} as primary key")

    def create_table(self, table_name: str) -> None:
        """
        Create table in database
        """
        engine = self.init_db_engine()
        with engine.connect() as con:
            con.execute(text(f"CREATE TABLE {table_name} ();"))
            con.commit()
        print(f"Table {table_name} created")

    def alter_table_add_column(self, table_name: str, column_name: str, data_type: str):
        """
        Add column to table
        """
        engine = self.init_db_engine()
        with engine.connect() as con:
            con.execute(text(f"ALTER TABLE {table_name} ADD COLUMN {column_name} {data_type};"))
            con.commit()
        print(f"Column {column_name} of type {data_type} added to table {table_name}")

    def alter_table_add_foreign_key(self, table_name: str, constraint_key: str, column_name: str,
                                    ref_table: str, ref_column: str) -> None:
        """
        Add foreign key to table column
        """
        engine = self.init_db_engine()
        with engine.connect() as con:
            con.execute(text(f"ALTER TABLE {table_name} ADD CONSTRAINT {constraint_key} "
                             f"FOREIGN KEY ({column_name}) REFERENCES {ref_table} ({ref_column});"))
            con.commit()
        print(f"Foreign key added to column {column_name} in table {table_name} referencing {ref_table}({ref_column})")

    def insert_data_into_compounds_targets(self):
        """
        Insert data into table compounds_targets
        """
        engine = self.init_db_engine()
        with engine.connect() as con:
            con.execute(text(f"INSERT INTO compounds_targets (compound_id, target_id) "
                             f"select compounds.compound_id, targets.target_id from tmp_comp_targets "
                             f"INNER JOIN compounds on tmp_comp_targets.smiles = compounds.smiles "
                             f"INNER JOIN targets on tmp_comp_targets.target = targets.sc_target_id"))
            con.commit()
        print(f"Data inserted into table compounds_targets")

    def insert_data_into_normalised_probabilities(self):
        """
        Insert data into table normalised_probabilities
        """
        engine = self.init_db_engine()
        with engine.connect() as con:
            con.execute(text(f"INSERT INTO normalised_probabilities (fingerprint_id, target_id, normalised_probability) "
                             f"select fingerprints.fingerprint_id, targets.target_id, tmp_target_normp.norm_prob "
                             f"from tmp_target_normp "
                             f"INNER JOIN fingerprints on tmp_target_normp.fingerprint = fingerprints.feature "
                             f"INNER JOIN targets on tmp_target_normp.target = targets.sc_target_id"))
            con.commit()
        print(f"Data inserted into table normalised_probabilities")

    def run_query(self, query):
        """
        Run query and return result
        """
        engine = self.init_db_engine()
        with engine.connect() as con:
            result = con.execute(text(query))
            return result.fetchall()


