// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamers/SingleLigandRotamerLibrary_and_StoredRotamerLibrarySpecification.util.hh
/// @brief  This file contains shared input strings and "right answer" code for the above classes.  (We are using strings instead of external files because we wish to test input from istreams (istringstreams) rather than from file objects.  (Admittedly this may be overzealous).
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_test_core_pack_rotamers_SingleLigandRotamerLibrary_and_StoredRotamerLibrarySpecification_hh
#define INCLUDED_test_core_pack_rotamers_SingleLigandRotamerLibrary_and_StoredRotamerLibrarySpecification_hh


// Test headers
#include <core/pack/rotamers/SingleLigandRotamerLibrary.hh> //for NamePosMap
#include <cxxtest/TestSuite.h>

//typedef std::map< std::string, core::Vector > NamePosMap;
//typedef numeric::xyzVector< core::Real > Vector; //probable
using namespace core::pack::rotamers; //for NamePosMap
typedef utility::vector1< NamePosMap > AtomPositions;

//utility functions for test data
//create the atom_positions comparator data structure.  Created straight from
//the pasted PDB; praise be to multiline editing
inline void generate_comparator( AtomPositions & return_vector ) {

	return_vector.clear();

	NamePosMap temp_map;

	temp_map.clear();
	temp_map[" C1 "] =  core::Vector(45.449, 15.367, -10.829);
	temp_map[" C2 "] =  core::Vector(46.057, 17.880, -11.049);
	temp_map[" C3 "] =  core::Vector(46.479, 16.487, -10.612);
	temp_map[" O1 "] =  core::Vector(45.033, 15.205, -12.009);
	temp_map[" O2 "] =  core::Vector(49.771, 14.920, -11.794);
	temp_map[" O3 "] =  core::Vector(47.856, 13.764, -10.551);
	temp_map[" O4 "] =  core::Vector(45.144, 14.729,  -9.779);
	temp_map[" O5 "] =  core::Vector(49.164, 15.658,  -9.421);
	temp_map[" O6 "] =  core::Vector(45.212, 18.455, -10.063);
	temp_map[" O7 "] =  core::Vector(47.639, 16.135, -11.364);
	temp_map[" P1 "] =  core::Vector(48.686, 15.032, -10.727);
	temp_map[" H1 "] =  core::Vector(45.505, 17.832, -11.995);
	temp_map[" H2 "] =  core::Vector(46.921, 18.534, -11.208);
	temp_map[" H3 "] =  core::Vector(46.747, 16.500,  -9.548);
	temp_map[" H4 "] =  core::Vector(44.968, 19.339, -10.383);
	return_vector.push_back(temp_map);

	temp_map.clear();
	temp_map[" C1 "] =  core::Vector(45.330, 14.633, -10.666);
	temp_map[" C2 "] =  core::Vector(45.937, 17.146, -10.887);
	temp_map[" C3 "] =  core::Vector(46.360, 15.753, -10.450);
	temp_map[" O1 "] =  core::Vector(44.914, 14.471, -11.847);
	temp_map[" O2 "] =  core::Vector(49.935, 15.492, -11.806);
	temp_map[" O3 "] =  core::Vector(49.210, 15.686,  -9.360);
	temp_map[" O4 "] =  core::Vector(45.025, 13.995,  -9.617);
	temp_map[" O5 "] =  core::Vector(48.728, 17.599, -10.998);
	temp_map[" O6 "] =  core::Vector(44.572, 17.356, -10.558);
	temp_map[" O7 "] =  core::Vector(47.520, 15.401, -11.202);
	temp_map[" P1 "] =  core::Vector(48.958, 16.103, -10.806);
	temp_map[" H1 "] =  core::Vector(46.055, 17.258, -11.971);
	temp_map[" H2 "] =  core::Vector(46.546, 17.921, -10.410);
	temp_map[" H3 "] =  core::Vector(46.628, 15.766,  -9.386);
	temp_map[" H4 "] =  core::Vector(44.346, 18.253, -10.854);
	return_vector.push_back(temp_map);

	temp_map.clear();
	temp_map[" C1 "] =  core::Vector(45.192, 14.638, -10.508);
	temp_map[" C2 "] =  core::Vector(45.800, 17.150, -10.729);
	temp_map[" C3 "] =  core::Vector(46.223, 15.758, -10.291);
	temp_map[" O1 "] =  core::Vector(45.651, 13.495, -10.782);
	temp_map[" O2 "] =  core::Vector(49.798, 15.497, -11.648);
	temp_map[" O3 "] =  core::Vector(49.073, 15.690,  -9.202);
	temp_map[" O4 "] =  core::Vector(43.989, 15.005, -10.375);
	temp_map[" O5 "] =  core::Vector(48.591, 17.604, -10.840);
	temp_map[" O6 "] =  core::Vector(45.970, 17.282, -12.132);
	temp_map[" O7 "] =  core::Vector(47.383, 15.406, -11.044);
	temp_map[" P1 "] =  core::Vector(48.821, 16.107, -10.648);
	temp_map[" H1 "] =  core::Vector(46.413, 17.910, -10.231);
	temp_map[" H2 "] =  core::Vector(44.755, 17.354, -10.472);
	temp_map[" H3 "] =  core::Vector(46.491, 15.771,  -9.228);
	temp_map[" H4 "] =  core::Vector(45.689, 18.182, -12.366);
	return_vector.push_back(temp_map);

}

//utility function: our string, hardcoded to generate an istringstream from.
inline std::string conformers_pdb_string() {
	return
		"HETATM    1  C1  4T1 X   1      45.449  15.367 -10.829  1.00 20.00           C  \n"
		"HETATM    2  C2  4T1 X   1      46.057  17.880 -11.049  1.00 20.00           C  \n"
		"HETATM    3  C3  4T1 X   1      46.479  16.487 -10.612  1.00 20.00           C  \n"
		"HETATM    4  O1  4T1 X   1      45.033  15.205 -12.009  1.00 20.00           O  \n"
		"HETATM    5  O2  4T1 X   1      49.771  14.920 -11.794  1.00 20.00           O  \n"
		"HETATM    6  O3  4T1 X   1      47.856  13.764 -10.551  1.00 20.00           O  \n"
		"HETATM    7  O4  4T1 X   1      45.144  14.729  -9.779  1.00 20.00           O  \n"
		"HETATM    8  O5  4T1 X   1      49.164  15.658  -9.421  1.00 20.00           O  \n"
		"HETATM    9  O6  4T1 X   1      45.212  18.455 -10.063  1.00 20.00           O  \n"
		"HETATM   10  O7  4T1 X   1      47.639  16.135 -11.364  1.00 20.00           O  \n"
		"HETATM   11  P1  4T1 X   1      48.686  15.032 -10.727  1.00 20.00           P  \n"
		"HETATM   12  H1  4T1 X   1      45.505  17.832 -11.995  1.00 20.00           H  \n"
		"HETATM   13  H2  4T1 X   1      46.921  18.534 -11.208  1.00 20.00           H  \n"
		"HETATM   14  H3  4T1 X   1      46.747  16.500  -9.548  1.00 20.00           H  \n"
		"HETATM   15  H4  4T1 X   1      44.968  19.339 -10.383  1.00 20.00           H  \n"
		"TER                                                                             \n"
		"HETATM    1  C1  4T1 X   1      45.330  14.633 -10.666  1.00 20.00           C  \n"
		"HETATM    2  C2  4T1 X   1      45.937  17.146 -10.887  1.00 20.00           C  \n"
		"HETATM    3  C3  4T1 X   1      46.360  15.753 -10.450  1.00 20.00           C  \n"
		"HETATM    4  O1  4T1 X   1      44.914  14.471 -11.847  1.00 20.00           O  \n"
		"HETATM    5  O2  4T1 X   1      49.935  15.492 -11.806  1.00 20.00           O  \n"
		"HETATM    6  O3  4T1 X   1      49.210  15.686  -9.360  1.00 20.00           O  \n"
		"HETATM    7  O4  4T1 X   1      45.025  13.995  -9.617  1.00 20.00           O  \n"
		"HETATM    8  O5  4T1 X   1      48.728  17.599 -10.998  1.00 20.00           O  \n"
		"HETATM    9  O6  4T1 X   1      44.572  17.356 -10.558  1.00 20.00           O  \n"
		"HETATM   10  O7  4T1 X   1      47.520  15.401 -11.202  1.00 20.00           O  \n"
		"HETATM   11  P1  4T1 X   1      48.958  16.103 -10.806  1.00 20.00           P  \n"
		"HETATM   12  H1  4T1 X   1      46.055  17.258 -11.971  1.00 20.00           H  \n"
		"HETATM   13  H2  4T1 X   1      46.546  17.921 -10.410  1.00 20.00           H  \n"
		"HETATM   14  H3  4T1 X   1      46.628  15.766  -9.386  1.00 20.00           H  \n"
		"HETATM   15  H4  4T1 X   1      44.346  18.253 -10.854  1.00 20.00           H  \n"
		"TER                                                                             \n"
		"HETATM    1  C1  4T1 X   1      45.192  14.638 -10.508  1.00 20.00           C  \n"
		"HETATM    2  C2  4T1 X   1      45.800  17.150 -10.729  1.00 20.00           C  \n"
		"HETATM    3  C3  4T1 X   1      46.223  15.758 -10.291  1.00 20.00           C  \n"
		"HETATM    4  O1  4T1 X   1      45.651  13.495 -10.782  1.00 20.00           O  \n"
		"HETATM    5  O2  4T1 X   1      49.798  15.497 -11.648  1.00 20.00           O  \n"
		"HETATM    6  O3  4T1 X   1      49.073  15.690  -9.202  1.00 20.00           O  \n"
		"HETATM    7  O4  4T1 X   1      43.989  15.005 -10.375  1.00 20.00           O  \n"
		"HETATM    8  O5  4T1 X   1      48.591  17.604 -10.840  1.00 20.00           O  \n"
		"HETATM    9  O6  4T1 X   1      45.970  17.282 -12.132  1.00 20.00           O  \n"
		"HETATM   10  O7  4T1 X   1      47.383  15.406 -11.044  1.00 20.00           O  \n"
		"HETATM   11  P1  4T1 X   1      48.821  16.107 -10.648  1.00 20.00           P  \n"
		"HETATM   12  H1  4T1 X   1      46.413  17.910 -10.231  1.00 20.00           H  \n"
		"HETATM   13  H2  4T1 X   1      44.755  17.354 -10.472  1.00 20.00           H  \n"
		"HETATM   14  H3  4T1 X   1      46.491  15.771  -9.228  1.00 20.00           H  \n"
		"HETATM   15  H4  4T1 X   1      45.689  18.182 -12.366  1.00 20.00           H  \n"
		"TER                                                                             \n"
		"REF_EN 77.0                                                                     \n"
		;
}

//note the stringified params does NOT have the PDB_ROTAMERS line
inline std::string params_string() {
	return
		"NAME 5T1\n" //tweaked name helps track that the code is using this string as desired
		"IO_STRING 4T1 1\n"
		"TYPE LIGAND\n"
		"AA UNK\n"
		"ATOM  C3  CH1   X   0.13\n"
		"ATOM  C1  COO   X   0.90\n"
		"ATOM  O1  OOC   X   -0.92\n"
		"ATOM  O4  OOC   X   -0.91\n"
		"ATOM  C2  CH2   X   0.11\n"
		"ATOM  O6  OH    X   -0.70\n"
		"ATOM  H4  Hpol  X   0.44\n"
		"ATOM  H1  Hapo  X   -0.05\n"
		"ATOM  H2  Hapo  X   0.06\n"
		"ATOM  O7  OH    X   -0.56\n"
		"ATOM  P1  Phos  X   1.38\n"
		"ATOM  O2  OOC   X   -0.95\n"
		"ATOM  O3  OOC   X   -0.95\n"
		"ATOM  O5  OOC   X   -1.01\n"
		"ATOM  H3  Hapo  X   0.02\n"
		"BOND  C1   C3 \n"
		"BOND  C1   O1 \n"
		"BOND  C1   O4 \n"
		"BOND  C2   C3 \n"
		"BOND  C2   O6 \n"
		"BOND  C3   O7 \n"
		"BOND  O2   P1 \n"
		"BOND  O3   P1 \n"
		"BOND  O5   P1 \n"
		"BOND  O7   P1 \n"
		"BOND  C2   H1 \n"
		"BOND  C2   H2 \n"
		"BOND  C3   H3 \n"
		"BOND  O6   H4 \n"
		"CHI 1  C3   C2   O6   H4 \n"
		"PROTON_CHI 1 SAMPLES 3 60 -60 180 EXTRA 1 20\n"
		"CHI 2  C2   C3   C1   O1 \n"
		"CHI 3  C1   C3   C2   O6 \n"
		"CHI 4  C1   C3   O7   P1 \n"
		"CHI 5  C3   O7   P1   O2 \n"
		"NBR_ATOM  C3 \n"
		"NBR_RADIUS 4.656937\n"
		"ICOOR_INTERNAL    C3     0.000000    0.000000    0.000000   C3    C1    O1 \n"
		"ICOOR_INTERNAL    C1     0.000000  180.000000    1.542636   C3    C1    O1 \n"
		"ICOOR_INTERNAL    O1     0.000000   60.982867    1.267998   C1    C3    O1 \n"
		"ICOOR_INTERNAL    O4  -174.390787   63.009576    1.260402   C1    C3    O1 \n"
		"ICOOR_INTERNAL    C2    91.227654   77.013363    1.521803   C3    C1    O1 \n"
		"ICOOR_INTERNAL    O6    -9.934447   69.199382    1.434911   C2    C3    C1 \n"
		"ICOOR_INTERNAL    H4   171.503849   70.501670    0.949990   O6    C2    C3 \n"
		"ICOOR_INTERNAL    H1   120.000660   70.971505    1.099989   C2    C3    O6 \n"
		"ICOOR_INTERNAL    H2   120.553167   71.264921    1.100021   C2    C3    H1 \n"
		"ICOOR_INTERNAL    O7  -124.405442   65.647991    1.456160   C3    C1    C2 \n"
		"ICOOR_INTERNAL    P1  -162.769241   57.458052    1.669434   O7    C3    C1 \n"
		"ICOOR_INTERNAL    O2  -138.955621   72.771120    1.528841   P1    O7    C3 \n"
		"ICOOR_INTERNAL    O3  -119.073941   71.472690    1.526110   P1    O7    O2 \n"
		"ICOOR_INTERNAL    O5  -117.943740   74.538247    1.531344   P1    O7    O3 \n"
		"ICOOR_INTERNAL    H3  -114.150353   67.608973    1.099991   C3    C1    O7 \n"
		;
}

//utility function: TS_ASSERT ALL OF THE THINGS!
//you can compare AtomPositions directly with TS_ASSERT, but I suspect it's a bad idea
inline void compare_AtomPositions( AtomPositions const & reference, AtomPositions const & data ){
	TS_ASSERT_EQUALS(reference.size(), data.size());

	//for each conformer
	for ( core::Size i(1); i<= reference.size(); ++i ) {
		NamePosMap const & reference_conformer(reference[i]);
		NamePosMap const & data_conformer     (data     [i]);

		TS_ASSERT_EQUALS(reference_conformer.size(), data_conformer.size());

		//extract the map keys (atom names) of the reference conformer
		//need the actual keys to co-iterate over both maps
		utility::vector1< std::string > keys;
		for ( auto const & it : reference_conformer ) {
			keys.push_back(it.first);
		}

		//get the data keys too to ensure they match
		utility::vector1< std::string > keys_data;
		for ( auto const & it : data_conformer ) {
			keys_data.push_back(it.first);
		}

		//ensure key sets match (ensure both NamePosMaps have the same atoms in them)
		TS_ASSERT_EQUALS(keys.size(), keys_data.size());
		for ( core::Size k(1); k<=keys.size(); ++k ) {
			TS_ASSERT_EQUALS(keys[k], keys_data[k]);
		}

		//for each atom in each conformer, ensure coordinates match
		for ( std::string const & key : keys ) {

			TS_ASSERT_DELTA(reference_conformer.at(key).x(), data_conformer.at(key).x(), 0.000001);
			TS_ASSERT_DELTA(reference_conformer.at(key).y(), data_conformer.at(key).y(), 0.000001);
			TS_ASSERT_DELTA(reference_conformer.at(key).z(), data_conformer.at(key).z(), 0.000001);

		}

	}
}

#endif //INCLUDED_test_core_pack_rotamers_SingleLigandRotamerLibrary_and_StoredRotamerLibrarySpecification_hh
