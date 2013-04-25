Generated HBond parameters by Matthew O'Meara in 2010


Adjusting hydrogen bond potentials for hydroxyl donors (newOH_params)
New OH potentials:
These are modeled off of the poly_AHdist_4 potential just shifted closer
poly_AHdist_4 has a peak at 1.902
poly_AHdist_dOH_1 has a peak at 1.735
poly_AHdist_dOH_2 has a peak at 1.753
poly_AHdist_dOH_3 has a peak at 1.84
AHX donors with PBA, HXL, CXA or CXL acceptors have polynomial poly_AHdist_dOH_1
HXL donors with PBA, AHX, HXL, CXA or CXL acceptors and dAHXaAHX have polynomial poly_AHdist_dOH_2
AHX or HXL donors to IMD or IME acceptors have polynomial poly_AHdist_dOH_3 

For each hydrogen bond class that uses the new OH parameters, shift
the fading functions in closer as well.



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



