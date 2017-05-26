2013/2/12
Parameter set "OLF_params_11c".

This parameter set has been tuned to work along with the hack-elec term.
If the sp2_params parameter set is used with hack-elec, the hydrogen bond
distributions are distorted, as would be expected when a new term is added
that dramatically alters the energy landscape in the region where hydrogen
bonds are formed.  To create these parameters, the sp2_params (9g) distance
polynomials were taken as a reference, and new polynomials were created so
that the total-energy derivatives in the region where hydrogen bonds
are formed matched; that is
d ( new hbond + elec ) / d dist = d ( old hbond ) / d dist

The paramter set has also been modified so that each of the polynomials range
from [-0.5,1.1] so that any one of the three polynomials can produce an
energy of +0.1.  The addition of the flag -corrections::score::hb_fade_energy
then smooths the hbond energy so that it smoothly turns off, eliminating
a previous derivative discontinuity.

This parameter set is intended for use with a fa_elec weight of 0.7
and the following set of flags on the command line.

-score::hb_sp2_chipen
-hbond_measure_sp3acc_BAH_from_hvy
-lj_hbond_hdis 1.75
-lj_hbond_OH_donor_dis 2.6
-hbond_params sp2_elec_params
-corrections::score::hb_fade_energy
-corrections::score::hb_sp2_outer_width .357
-corrections::chemical::expand_st_chi2sampling
-smooth_fa_elec
-elec_min_dis 1.6
-elec_r_option false


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


