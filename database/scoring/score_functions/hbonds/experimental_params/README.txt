Generated HBond parameters by Matthew O'Meara in 2011

This parameter set contains the following corrections:
  - helix_hb_06_2009 correction
  - newOH correction
  - BBasSC correction

NOTE: It is expected that when using the newOH correction the following flags are also used:
  -lj_hbond_hdis 1.75
  -lj_hbond_OH_donor_dis 2.6


README.txt: this file
schema.sql: sql schema definition for tables
HBWeightType.csv: Define the components of the overall score function where hydrogen bond scores will be recorded
HBDonChemType.csv: Define chemical classes of donors
HBAccChemType.csv: Define chemical classes of acceptors
HBAccHybridization.csv: Geometric features cean be defined in terms of the hybridization of the acceptor moity.
HybridizationType.csv: Define chemical classes for orbital hybridization
HBSeqSep.csv: An identification feature of an hbond (proxy for secondary structure) 
HBPoly1D.csv: 1d Polynomials for the Hydrogen Bond Score Term
HBEval.csv: Evaluation parameters for each specific hydrogen bond type
HBFadeIntervals.csv fade intervals used simplified cross terms in the hbond evaluation.

To create the database execute

sqlite3 minirosetta_database.db3 < schema.sql

Reading in the schema file will load all csv data files into the the database.


# Utility scripts
sqlite3_io.py: script to transform flatfile database data to sqlite database data and visa-versa 
convert_names_to_ids.py



