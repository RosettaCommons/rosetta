Generated HBond parameters by Matthew O'Meara in 2010

This is a correction made by Phil to the hbonds with ring acceptors:

///////////////////////////////////////////

 /**
  /// this is the old polynomial: which goes down to -0.22
  create_poly7(POLY_xHRing, MIN_xH, MAX_xH, // 0.7608, 1.089,e
  37.744316,-117.731674,143.0759275,-86.2258835,26.7448175,-4.4699705,0.6458455)

  /// this was my original fix to make it go to -0.5:
  create_poly7(POLY_xHRing, MIN_xH, MAX_xH, // 0.7608, 1.089,e
  37.744316,-117.731674,143.0759275,-86.2258835,26.7448175,-4.4699705,0.3658455)
 **/
 /// Phil:
 /// this is a new fix, which goes to -0.5 while preserving the 0-crossing point of the old polynomial (~140 degrees)
 /// the 0-crossing with the previous fix was ~108 degrees
 /// this makes it much more restrictive than the previous fix
 /// I'm making this change in response to bad geometry at the ring acceptor of protein-DNA interfaces
 ///

static Real const xHRing_hack_factor( 0.5 / 0.22 );
create_poly7(POLY_xHRing_Phil, MIN_xH, MAX_xH, // 0.7608, 1.089,e
       	     xHRing_hack_factor * 37.744316, xHRing_hack_factor * -117.731674, xHRing_hack_factor * 143.0759275,
	     xHRing_hack_factor * -86.2258835, xHRing_hack_factor * 26.7448175, xHRing_hack_factor * -4.4699705,
	     xHRing_hack_factor * 0.6458455 )

In the old evaluation case switch there was:

  case hbe_RINGSC:
    POLY_AHdisSP2(dAHdis, Pr,dPr);
    POLY_xDSP2short(dxD, PSxD,dPSxD);
    if (core::options::option[core::options::OptionKeys::corrections::score::hbond_His_Phil_fix]) {
      POLY_xHRing_Phil(dxH, PSxH,dPSxH);
    } else {
      POLY_xHRing(dxH, PSxH,dPSxH);
    }
    POLY_xDSP2long(dxD, PLxD,dPLxD);
    PLxH = PSxH; dPLxH = dPSxH;
    break;
  }


//////////////////////////////////////


now to use this correction put 

     -scoring:hbond_params His_Phil_fix 

on the command line to use the parameters in this directory.







Directory Contents:

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



