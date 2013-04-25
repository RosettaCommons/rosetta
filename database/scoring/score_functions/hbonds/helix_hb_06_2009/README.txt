Generated HBond parameters by Matthew O'Meara in 2010

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





