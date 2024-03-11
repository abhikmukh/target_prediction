#!/usr/bin/env python

import psycopg2

conn = psycopg2.connect("host=localhost dbname=nb_target user=postgres password=abhik1234")


cur = conn.cursor()


# create compounds table
cur.execute("""CREATE TABLE compounds(
    monomer_id text,
    smiles text
)
""")

conn.commit()


with open('map_compounds.csv', 'r') as f:
    # Notice that we don't need the `csv` module.
    next(f) # Skip the header row.
    cur.copy_from(f, 'compounds', sep=',')

conn.commit()

cur.execute("""ALTER TABLE compounds ADD COLUMN compound_id BIGSERIAL PRIMARY KEY""")
conn.commit()
print("finish creating compounds table")

# create fingerprint table
cur.execute("""CREATE TABLE fingerprints(
    feature bigint,
    frequency_all bigint
)
""")
conn.commit()


with open('map_fingerprints.csv', 'r') as f:
    # Notice that we don't need the `csv` module.
    cur.copy_from(f, 'fingerprints', sep=',')

conn.commit()

cur.execute("""ALTER TABLE fingerprints ADD COLUMN fingerprint_id BIGSERIAL PRIMARY KEY""")
conn.commit()

print("finish creating fingerprint table")

# create patents table
cur.execute("""CREATE TABLE targets(
    sc_target_id text,
    laplaciank numeric
)
""")

conn.commit()

with open('map_target_laplacian.csv', 'r') as f:
    next(f)
    # Notice that we don't need the `csv` module.
    cur.copy_from(f, 'targets', sep=',')

conn.commit()

cur.execute("""ALTER TABLE targets ADD COLUMN target_id BIGSERIAL PRIMARY KEY""")

conn.commit()
print("finish creating targets table")


# create temporary table tmp_comp_patents
cur.execute("""CREATE TABLE tmp_comp_targets (
  smiles text,
  target text
)""")


conn.commit()

with open('map_compounds_targets.csv', 'r') as f:
    next(f)
    # Notice that we don't need the `csv` module.
    cur.copy_from(f, 'tmp_comp_targets', sep=',')

conn.commit()
print("finish creating tmp_comp table")


# create empty table compounds_patents
cur.execute("""CREATE TABLE compounds_targets ()""")


conn.commit()


cur.execute("""alter table compounds_targets add column compound_id bigint""")


conn.commit()


cur.execute("""alter table compounds_targets add column target_id bigint""")

conn.commit()


cur.execute("""alter table compounds_targets
                add constraint fk_compounds_compounds_targets foreign key (compound_id)
                references compounds (compound_id)
            """)


conn.commit()


cur.execute("""alter table compounds_targets
                add constraint fk_targets_compounds_targets foreign key (target_id)
                references targets (target_id)
            """)

conn.commit()


cur.execute("""INSERT INTO compounds_targets (compound_id, target_id) 
                select compounds.compound_id, targets.target_id from tmp_comp_targets 
INNER JOIN compounds on tmp_comp_targets.smiles = compounds.smiles
INNER JOIN targets on tmp_comp_targets.target = targets.sc_target_id
                """)


conn.commit()

# create tempoprary table tmp_target_normp
cur.execute("""CREATE TABLE tmp_target_normp(
    target text,
    fingerprint bigint,
    norm_prob numeric
)
""")

conn.commit()

with open('map_target_norm_p.csv', 'r') as f:
    next(f)
    # Notice that we don't need the `csv` module.
    cur.copy_from(f, 'tmp_target_normp', sep=',')

conn.commit()

cur.execute("""CREATE TABLE normalised_probablities (
  target_id bigint,
  fingerprint_id bigint
)""")

conn.commit()

cur.execute("""alter table normalised_probablities
                add constraint fk_fingerprints_normalised_probablities foreign key (fingerprint_id)
                references fingerprints (fingerprint_id)
            """)


conn.commit()


cur.execute("""alter table normalised_probablities
                add constraint fk_targets_normalised_probablities foreign key (target_id)
                references targets (target_id)
            """)


conn.commit()

cur.execute("""alter table normalised_probablities add column normalised_probablity numeric""")


conn.commit()


cur.execute("""INSERT INTO normalised_probablities (fingerprint_id, target_id, normalised_probablity)
                select fingerprints.fingerprint_id, targets.target_id, tmp_target_normp.norm_prob from tmp_target_normp 
                INNER JOIN fingerprints on tmp_target_normp.fingerprint = fingerprints.feature
                INNER JOIN targets on tmp_target_normp.target = targets.sc_target_id""")


conn.commit()


