import os

from .database_utils import DatabaseConnector


def run(data_dir: str = "data") -> None:
    base_path = os.getcwd()
    print(base_path)
    database_connector = DatabaseConnector(base_path=base_path, cred_file="db_credential.yaml")
    database_connector.upload_to_db(csv_file=(os.path.join(data_dir, "map_compounds.csv")), table_name="compounds")
    database_connector.alter_table_data_type(table_name="compounds", column_name="monomer_id", data_type="TEXT")
    database_connector.alter_table_data_type(table_name="compounds", column_name="smiles", data_type="TEXT")
    database_connector.add_column_with_primary_key(table_name="compounds", column_name="compound_id",
                                                   data_type="BIGSERIAL")
    print(f"Data uploaded to table compounds")

    database_connector.upload_to_db(csv_file=(os.path.join(data_dir, "map_fingerprints.csv")),
                                    table_name="fingerprints")
    database_connector.alter_table_data_type(table_name="fingerprints", column_name="feature", data_type="BIGINT")
    database_connector.alter_table_data_type(table_name="fingerprints", column_name="frequency_all", data_type="BIGINT")
    database_connector.add_column_with_primary_key(table_name="fingerprints", column_name="fingerprint_id",
                                                   data_type="BIGSERIAL")
    print(f"Data uploaded to table fingerprints")

    database_connector.upload_to_db(csv_file=(os.path.join(data_dir, "map_target_laplacian.csv")), table_name="targets")
    database_connector.add_column_with_primary_key(table_name="targets", column_name="target_id", data_type="BIGSERIAL")
    print(f"Data uploaded to table targets")

    database_connector.upload_to_db(csv_file=(os.path.join(data_dir, "map_compounds_targets.csv")),
                                    table_name="tmp_comp_targets")
    database_connector.alter_table_data_type(table_name="tmp_comp_targets", column_name="smiles", data_type="TEXT")
    database_connector.alter_table_data_type(table_name="tmp_comp_targets", column_name="target", data_type="TEXT")
    print(f"Data uploaded to table tmp_comp_targets")

    database_connector.create_table(table_name="compounds_targets")
    database_connector.alter_table_add_column(table_name="compounds_targets", column_name="compound_id",
                                              data_type="BIGINT")
    database_connector.alter_table_add_column(table_name="compounds_targets", column_name="target_id",
                                              data_type="BIGINT")
    database_connector.alter_table_add_foreign_key(table_name="compounds_targets",
                                                   constraint_key="fk_compounds_compounds_targets",
                                                   column_name="compound_id", ref_table="compounds",
                                                   ref_column="compound_id")
    database_connector.alter_table_add_foreign_key(table_name="compounds_targets",
                                                   constraint_key="fk_targets_compounds_targets",
                                                   column_name="target_id", ref_table="targets",
                                                   ref_column="target_id")
    database_connector.insert_data_into_compounds_targets()

    print(f"Data uploaded to table compounds_targets")

    database_connector.upload_to_db(csv_file=(os.path.join(data_dir, "map_target_norm_p.csv")),
                                    table_name="tmp_target_normp")
    database_connector.alter_table_data_type(table_name="tmp_target_normp", column_name="target", data_type="TEXT")
    database_connector.alter_table_data_type(table_name="tmp_target_normp", column_name="fingerprint",
                                             data_type="BIGINT")
    database_connector.alter_table_data_type(table_name="tmp_target_normp", column_name="norm_prob",
                                             data_type="NUMERIC")
    print(f"Data uploaded to table tmp_target_normp")

    database_connector.create_table(table_name="normalised_probabilities")
    database_connector.alter_table_add_column(table_name="normalised_probabilities", column_name="target_id",
                                              data_type="BIGINT")
    database_connector.alter_table_add_column(table_name="normalised_probabilities", column_name="fingerprint_id",
                                              data_type="BIGINT")
    database_connector.alter_table_add_foreign_key(table_name="normalised_probabilities",
                                                   constraint_key="fk_fingerprints_normalised_probabilities",
                                                   column_name="fingerprint_id", ref_table="fingerprints",
                                                   ref_column="fingerprint_id")
    database_connector.alter_table_add_foreign_key(table_name="normalised_probabilities",
                                                   constraint_key="fk_targets_normalised_probabilities",
                                                   column_name="target_id", ref_table="targets", ref_column="target_id")
    database_connector.alter_table_add_column(table_name="normalised_probabilities",
                                              column_name="normalised_probability", data_type="NUMERIC")
    database_connector.insert_data_into_normalised_probabilities()
    print(f"Data uploaded to table normalised_probabilities")


if __name__ == "__main__":
    run()
