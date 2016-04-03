// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/core/import_pose/3-letter_code_ambiguity_resolution.cxxtest.hh
/// @brief   Test suite for ensuring that a .pdb file with shared or ambiguous 3-letter codes can be built into a Pose.
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package header
#include <core/import_pose/import_pose.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// C++ header
#include <string>


class ThreeLetterCodeAmbiguityTests : public CxxTest::TestSuite {
public: // Standard methods ///////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init();
	}

	// Destruction
	void tearDown()
	{}


public: // Tests //////////////////////////////////////////////////////////////
	// Confirm that two different ResidueTypes are created from .pdb file with
	// only a single 3-letter code used.
	void test_resolution_of_shared_3_letter_codes()
	{
		using namespace core;

		std::string const flavins(
			"HETNAM     Fl  X   1  flavin\n"
			"HETNAM     Fl  Y   1  flavin_dihydride\n"
			"HETATM    1  N5  Fl  X   1       0.848  -0.156  -0.082  1.00  0.00           N  \n"
			"HETATM    2  C4a Fl  X   1       0.575  -1.348  -0.448  1.00  0.00           C  \n"
			"HETATM    3  C4  Fl  X   1       1.611  -2.357  -0.490  1.00  0.00           C  \n"
			"HETATM    4  O4  Fl  X   1       2.766  -2.064  -0.170  1.00  0.00           O  \n"
			"HETATM    5  N3  Fl  X   1       1.300  -3.593  -0.876  1.00  0.00           N  \n"
			"HETATM    6  C2  Fl  X   1      -0.029  -3.908  -1.238  1.00  0.00           C  \n"
			"HETATM    7  O2  Fl  X   1      -0.318  -5.063  -1.599  1.00  0.00           O  \n"
			"HETATM    8  N1  Fl  X   1      -1.040  -2.910  -1.194  1.00  0.00           N  \n"
			"HETATM    9  C10 Fl  X   1      -0.775  -1.717  -0.829  1.00  0.00           C  \n"
			"HETATM   10  N10 Fl  X   1      -1.770  -0.754  -0.791  1.00  0.00           N  \n"
			"HETATM   11  C9a Fl  X   1      -1.454   0.528  -0.392  1.00  0.00           C  \n"
			"HETATM   12  C9  Fl  X   1      -2.438   1.508  -0.347  1.00  0.00           C  \n"
			"HETATM   13  C8  Fl  X   1      -2.111   2.797   0.056  1.00  0.00           C  \n"
			"HETATM   14  C7  Fl  X   1      -0.801   3.102   0.412  1.00  0.00           C  \n"
			"HETATM   15  C6  Fl  X   1       0.181   2.121   0.366  1.00  0.00           C  \n"
			"HETATM   16  C5a Fl  X   1      -0.145   0.832  -0.036  1.00  0.00           C  \n"
			"HETATM   17  C7a Fl  X   1      -0.441   4.489   0.847  1.00  0.00           C  \n"
			"HETATM   18  C8a Fl  X   1      -3.175   3.852   0.104  1.00  0.00           C  \n"
			"TER      19      Fl  X   1 \n"
			"HETATM   20  N5  Fl  Y   1       4.745 -12.417  -2.665  1.00  0.00           N  \n"
			"HETATM   21  C4a Fl  Y   1       4.420 -13.731  -3.074  1.00  0.00           C  \n"
			"HETATM   22  C4  Fl  Y   1       5.460 -14.753  -3.118  1.00  0.00           C  \n"
			"HETATM   23  O4  Fl  Y   1       6.613 -14.452  -2.796  1.00  0.00           O  \n"
			"HETATM   24  N3  Fl  Y   1       5.170 -16.000  -3.502  1.00  0.00           N  \n"
			"HETATM   25  C2  Fl  Y   1       3.858 -16.340  -3.868  1.00  0.00           C  \n"
			"HETATM   26  O2  Fl  Y   1       3.582 -17.501  -4.227  1.00  0.00           O  \n"
			"HETATM   27  N1  Fl  Y   1       2.855 -15.352  -3.825  1.00  0.00           N  \n"
			"HETATM   28  C10 Fl  Y   1       3.170 -14.041  -3.419  1.00  0.00           C  \n"
			"HETATM   29  N10 Fl  Y   1       2.169 -13.066  -3.379  1.00  0.00           N  \n"
			"HETATM   30  C9a Fl  Y   1       2.464 -11.771  -2.981  1.00  0.00           C  \n"
			"HETATM   31  C9  Fl  Y   1       1.465 -10.806  -2.943  1.00  0.00           C  \n"
			"HETATM   32  C8  Fl  Y   1       1.772  -9.511  -2.543  1.00  0.00           C  \n"
			"HETATM   33  C7  Fl  Y   1       3.076  -9.184  -2.182  1.00  0.00           C  \n"
			"HETATM   34  C6  Fl  Y   1       4.074 -10.151  -2.221  1.00  0.00           C  \n"
			"HETATM   35  C5a Fl  Y   1       3.763 -11.444  -2.621  1.00  0.00           C  \n"
			"HETATM   36  C7a Fl  Y   1       3.414  -7.791  -1.750  1.00  0.00           C  \n"
			"HETATM   37  C8a Fl  Y   1       0.693  -8.472  -2.502  1.00  0.00           C  \n"
			"TER      38      Fl  Y   1 \n"
			"END\n" );

		pose::Pose pose;
		import_pose::pose_from_pdbstring( pose, flavins );

		TS_ASSERT_EQUALS( pose.n_residue(), 2 );
		TS_ASSERT_EQUALS( pose.residue( 1 ).name(), "flavin" );
		TS_ASSERT_EQUALS( pose.residue( 2 ).name(), "flavin_dihydride" );
		TS_ASSERT( pose.residue( 1 ).natoms() != pose.residue( 2 ).natoms() );
	}
};  // class ThreeLetterCodeAmbiguityTests
