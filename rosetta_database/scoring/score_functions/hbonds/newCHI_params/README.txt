Generated HBond parameters by Matthew O'Meara and Andrew Leaver-Fay in 2011

Hydrogen bond parameter set that extends the cosBAH functions down to
-.3766.  Native structures have many examples of hydrogen bonds with
BAH angle > 90 (where 0 is linear). To recognize these hydrogen bonds
te the fade functions and the polynomials are extended. Note to use
this, currently requires adjusting the MIN_xH in
src/core/scoring/hbonds/constants.hh (7/12/11).



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



