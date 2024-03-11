from naive_bayes import create_training_and_test_set, generate_fp_data, generate_normalised_probabilities, \
    database_loader


def main():
    data_dir = "naive_bayes/data"
    train_test_prep = create_training_and_test_set.TrainTestPreparation(data_dir)
    train_test_prep.create_training_and_test_set("processed_smiles_target.csv")
    fp_generator = generate_fp_data.FingerprintGenerator(data_dir)
    fp_generator.generate_fingerprint_set("train.csv", "train_fp.csv")
    fp_generator.generate_fingerprint_set_new("train.csv", "target_fp.csv")
    generate_normalised_probabilities.run(data_dir=data_dir)
    database_loader.run(data_dir=data_dir)


if __name__ == "__main__":
    main()

