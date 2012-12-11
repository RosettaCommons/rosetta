Generated HBond parameters by Matthew O'Meara and Andrew Leaver-Fay in 2011-2012 

12/11/12: This is the 9g parameter set.  The 9a parameter set had problems with it resulting from
the conversion from the "multiplicative" form the sp2 potential to the "additive" form.

This parameter set is intended for use with the following additional flags:

-score::hb_sp2_chipen
-hb_sp2_BAH180_rise 0.75  <--- currently the default parameter for this option, so it may be omitted.
-hbond_measure_sp3acc_BAH_from_hvy
-lj_hbond_hdis 1.75
-lj_hbond_OH_donor_dis 2.6
-score:weights sp2_correction.wts

The sp2_correction.wts parameter set differs from score12 only in the addition of the "yhh_planarity"
term, which puts a sinusoidal penalty on tyrosine's chi2 when it deviates from 0 or 180, and
it also increases the weight on alpha-helical hydrogen bonds (hbond_sr_bb) so that they are weighted
equally to the beta-sheet hydrogen bonds (hbond_lr_bb).  Or rather, the new weight set does not
down-weight alpha-helical hydrogne bonds the way that score12 does. Alpha-helical hydrogen bonds
don't actually get stronger with the new weight set; the sp2 potential penalizes them significantly
in comparison to the closer-to-ideal beta-sheet hydrogen bonds -- the downweighting that score12
applies no longer is necessary.

...

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

