// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/import_pose/3-letter_code_ambiguity_resolution.cxxtest.hh
/// @brief   Test suite for ensuring that a .pdb file with shared or ambiguous 3-letter codes can be built into a Pose.
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package header
#include <core/import_pose/import_pose.hh>

// Project headers
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

// C++ header
#include <string>

static THREAD_LOCAL basic::Tracer TR("core.import_pose.3-letter_code_ambiguity_resolution.cxxtest");

class ThreeLetterCodeAmbiguityTests : public CxxTest::TestSuite {
public: // Standard methods ///////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init_with_additional_options( "-include_sugars -alternate_3_letter_codes pdb_sugar" );
	}

	// Destruction
	void tearDown()
	{}


public: // Tests //////////////////////////////////////////////////////////////
	// Confirm that two different ResidueTypes are created from .pdb file with
	// only a single 3-letter code used.
	void test_resolution_of_shared_3_letter_codes()
	{
#ifdef MULTI_THREADED
		try {
#else
		{
#endif
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

			TS_ASSERT_EQUALS( pose.size(), 2 );
			TS_ASSERT_EQUALS( pose.residue( 1 ).name(), "flavin" );
			TS_ASSERT_EQUALS( pose.residue( 2 ).name(), "flavin_dihydride" );
			TS_ASSERT( pose.residue( 1 ).natoms() != pose.residue( 2 ).natoms() );
#ifdef MULTI_THREADED
		} catch(utility::excn::EXCN_Base& excn) {
			std::string expected( "ERROR: Error in ScoringManager: the carbohydrate CHIEnergyFunction is fundamentally not threadsafe, and cannot be used in a multithreaded environment.  Please contact Jason Labonte (JWLabonte@jhu.edu) to complain about this." );
			TS_ASSERT_EQUALS( excn.msg().substr( excn.msg().find( "ERROR: " ), expected.size() ), expected );
		}
#else
		}
#endif
	}

	// Confirm that ResidueTypes are created from .pdb file with the proper
	// branching variants.
	void test_resolution_of_linkage_information()
	{
		using namespace core;
		using namespace conformation;

#ifdef MULTI_THREADED
		try {
#else
		{
#endif

			std::string const Lex_Rosetta_format(
				"HETNAM     Glc A   1  ->4)-beta-D-Glcp, 2-acetylamino-d-deoxy-\n"
				"HETNAM     Gal A   2  ->4)-beta-D-Galp\n"
				"HETNAM     Fuc B   1  ->4)-alpha-L-Fucp\n"
				"LINK         O3  Glc A   1                 C1  Fuc B   1     1555   1555  1.5   \n"
				"HETATM    1  C1  Glc A   1      35.710 122.693  36.907  1.00 43.08           C  \n"
				"HETATM    2  C2  Glc A   1      36.264 121.664  35.887  1.00 40.87           C  \n"
				"HETATM    3  C3  Glc A   1      35.396 120.362  35.896  1.00 39.23           C  \n"
				"HETATM    4  C4  Glc A   1      35.215 119.833  37.360  1.00 34.96           C  \n"
				"HETATM    5  C5  Glc A   1      34.756 121.002  38.316  1.00 37.04           C  \n"
				"HETATM    6  C6  Glc A   1      34.689 120.700  39.820  1.00 38.17           C  \n"
				"HETATM    7  CN2 Glc A   1      37.363 122.3    33.878  1.00 40.62           C  \n"
				"HETATM    8 CAN2 Glc A   1      37.043 123.112  32.473  1.00 38.75           C  \n"
				"HETATM    9  N2  Glc A   1      36.245 122.328  34.583  1.00 40.83           N  \n"
				"HETATM   10  O1  Glc A   1      36.369 123.964  36.940  1.00 49.71           O  \n"
				"HETATM   11  O3  Glc A   1      36.026 119.359  35.094  1.00 38.07           O  \n"
				"HETATM   12  O4  Glc A   1      34.235 118.792  37.370  1.00 31.94           O  \n"
				"HETATM   13  O5  Glc A   1      35.660 122.104  38.210  1.00 39.15           O  \n"
				"HETATM   14  O6  Glc A   1      35.928 120.224  40.287  1.00 40.77           O  \n"
				"HETATM   15 OCN2 Glc A   1      38.517 122.472  34.271  1.00 40.92           O  \n"
				"HETATM   16  C1  Gal A   2      34.615 117.665  38.155  1.00 30.24           C  \n"
				"HETATM   17  C2  Gal A   2      33.365 116.921  38.692  1.00 30.24           C  \n"
				"HETATM   18  C3  Gal A   2      33.836 115.668  39.463  1.00 27.25           C  \n"
				"HETATM   19  C4  Gal A   2      34.765 114.785  38.592  1.00 26.13           C  \n"
				"HETATM   20  C5  Gal A   2      35.931 115.629  38.011  1.00 26.72           C  \n"
				"HETATM   21  C6  Gal A   2      36.818 114.855  37.032  1.00 25.57           C  \n"
				"HETATM   22  O2  Gal A   2      32.626 117.740  39.575  1.00 34.02           O  \n"
				"HETATM   23  O3  Gal A   2      32.706 114.881  39.893  1.00 28.50           O  \n"
				"HETATM   24  O4  Gal A   2      33.936 114.234  37.591  1.00 27.66           O  \n"
				"HETATM   25  O5  Gal A   2      35.406 116.786  37.348  1.00 27.55           O  \n"
				"HETATM   26  O6  Gal A   2      37.611 113.891  37.692  1.00 25.73           O  \n"
				"HETATM   27  C1  Fuc B   1      35.546 119.238  33.750  1.00 36.63           C  \n"
				"HETATM   28  C2  Fuc B   1      36.615 118.464  32.939  1.00 32.99           C  \n"
				"HETATM   29  C3  Fuc B   1      36.681 116.991  33.421  1.00 30.04           C  \n"
				"HETATM   30  C4  Fuc B   1      35.284 116.331  33.366  1.00 23.12           C  \n"
				"HETATM   31  C5  Fuc B   1      34.240 117.183  34.142  1.00 29.38           C  \n"
				"HETATM   32  C6  Fuc B   1      32.809 116.640  34.044  1.00 27.79           C  \n"
				"HETATM   33  O2  Fuc B   1      37.839 119.137  33.073  1.00 35.20           O  \n"
				"HETATM   34  O3  Fuc B   1      37.527 116.202  32.617  1.00 23.08           O  \n"
				"HETATM   35  O4  Fuc B   1      34.915 116.247  31.999  1.00 20.40           O  \n"
				"HETATM   36  O5  Fuc B   1      34.280 118.546  33.673  1.00 31.96           O  \n"
				"TER      37      Fuc B   1\n"
				"END\n" );
			std::string const Lex_pdb_format(
				"HETNAM     NAG N-ACETYL-D-GLUCOSAMINE\n"
				"HETNAM     GAL BETA-D-GALACTOSE\n"
				"HETNAM     FUC ALPHA-L-FUCOSE\n"
				"LINK         O4  NAG A   1                 C1  GAL A   2     1555   1555  1.5   \n"
				"LINK         O3  NAG A   1                 C1  FUC A   3     1555   1555  1.5   \n"
				"HETATM    1  C1  NAG A   1      35.710 122.693  36.907  1.00 43.08           C  \n"
				"HETATM    2  C2  NAG A   1      36.264 121.664  35.887  1.00 40.87           C  \n"
				"HETATM    3  C3  NAG A   1      35.396 120.362  35.896  1.00 39.23           C  \n"
				"HETATM    4  C4  NAG A   1      35.215 119.833  37.360  1.00 34.96           C  \n"
				"HETATM    5  C5  NAG A   1      34.756 121.002  38.316  1.00 37.04           C  \n"
				"HETATM    6  C6  NAG A   1      34.689 120.700  39.820  1.00 38.17           C  \n"
				"HETATM    7  C7  NAG A   1      37.363 122.3    33.878  1.00 40.62           C  \n"
				"HETATM    8  C8  NAG A   1      37.043 123.112  32.473  1.00 38.75           C  \n"
				"HETATM    9  N2  NAG A   1      36.245 122.328  34.583  1.00 40.83           N  \n"
				"HETATM   10  O1  NAG A   1      36.369 123.964  36.940  1.00 49.71           O  \n"
				"HETATM   11  O3  NAG A   1      36.026 119.359  35.094  1.00 38.07           O  \n"
				"HETATM   12  O4  NAG A   1      34.235 118.792  37.370  1.00 31.94           O  \n"
				"HETATM   13  O5  NAG A   1      35.660 122.104  38.210  1.00 39.15           O  \n"
				"HETATM   14  O6  NAG A   1      35.928 120.224  40.287  1.00 40.77           O  \n"
				"HETATM   15  O7  NAG A   1      38.517 122.472  34.271  1.00 40.92           O  \n"
				"HETATM   16  C1  GAL A   2      34.615 117.665  38.155  1.00 30.24           C  \n"
				"HETATM   17  C2  GAL A   2      33.365 116.921  38.692  1.00 30.24           C  \n"
				"HETATM   18  C3  GAL A   2      33.836 115.668  39.463  1.00 27.25           C  \n"
				"HETATM   19  C4  GAL A   2      34.765 114.785  38.592  1.00 26.13           C  \n"
				"HETATM   20  C5  GAL A   2      35.931 115.629  38.011  1.00 26.72           C  \n"
				"HETATM   21  C6  GAL A   2      36.818 114.855  37.032  1.00 25.57           C  \n"
				"HETATM   22  O2  GAL A   2      32.626 117.740  39.575  1.00 34.02           O  \n"
				"HETATM   23  O3  GAL A   2      32.706 114.881  39.893  1.00 28.50           O  \n"
				"HETATM   24  O4  GAL A   2      33.936 114.234  37.591  1.00 27.66           O  \n"
				"HETATM   25  O5  GAL A   2      35.406 116.786  37.348  1.00 27.55           O  \n"
				"HETATM   26  O6  GAL A   2      37.611 113.891  37.692  1.00 25.73           O  \n"
				"HETATM   27  C1  FUC A   3      35.546 119.238  33.750  1.00 36.63           C  \n"
				"HETATM   28  C2  FUC A   3      36.615 118.464  32.939  1.00 32.99           C  \n"
				"HETATM   29  C3  FUC A   3      36.681 116.991  33.421  1.00 30.04           C  \n"
				"HETATM   30  C4  FUC A   3      35.284 116.331  33.366  1.00 23.12           C  \n"
				"HETATM   31  C5  FUC A   3      34.240 117.183  34.142  1.00 29.38           C  \n"
				"HETATM   32  C6  FUC A   3      32.809 116.640  34.044  1.00 27.79           C  \n"
				"HETATM   33  O2  FUC A   3      37.839 119.137  33.073  1.00 35.20           O  \n"
				"HETATM   34  O3  FUC A   3      37.527 116.202  32.617  1.00 23.08           O  \n"
				"HETATM   35  O4  FUC A   3      34.915 116.247  31.999  1.00 20.40           O  \n"
				"HETATM   36  O5  FUC A   3      34.280 118.546  33.673  1.00 31.96           O  \n"
				"TER      37      FUC A   3\n"
				"END\n" );

			pose::Pose poseA, poseB;

			TR << "Making Rosetta-format pose " << std::endl;
			import_pose::pose_from_pdbstring( poseA, Lex_Rosetta_format );
			Residue const & resA1( poseA.residue( 1 ) );
			Residue const & resA2( poseA.residue( 2 ) );
			Residue const & resA3( poseA.residue( 3 ) );
			TR << "Rosetta format: " << resA1.name() << " " << resA2.name() << " " << resA3.name() << std::endl;
			TS_ASSERT( resA1.is_lower_terminus() );
			TS_ASSERT( resA1.is_branch_point() );
			TS_ASSERT( resA2.is_upper_terminus() );
			TS_ASSERT( resA3.is_upper_terminus() );

			TR << "Making PDB-format pose " << std::endl;
			import_pose::pose_from_pdbstring( poseB, Lex_pdb_format );
			Residue const & resB1( poseB.residue( 1 ) );
			Residue const & resB2( poseB.residue( 2 ) );
			Residue const & resB3( poseB.residue( 3 ) );
			TR << "PDB format: " << resB1.name() << " " << resB2.name() << " " << resB3.name() << std::endl;
			TS_ASSERT( resB1.is_lower_terminus() );
			TS_ASSERT( resB1.is_branch_point() );
			TS_ASSERT( resB2.is_upper_terminus() );
			TS_ASSERT( resB3.is_upper_terminus() );

			TS_ASSERT_EQUALS( poseA.size(), poseB.size() );
			TS_ASSERT_EQUALS( poseA.fold_tree().to_string(), poseB.fold_tree().to_string() );

#ifdef MULTI_THREADED
		} catch(utility::excn::EXCN_Base& excn) {
			std::string expected( "ERROR: Error in ScoringManager: the carbohydrate CHIEnergyFunction is fundamentally not threadsafe, and cannot be used in a multithreaded environment.  Please contact Jason Labonte (JWLabonte@jhu.edu) to complain about this." );
			TS_ASSERT_EQUALS( excn.msg().substr( excn.msg().find( "ERROR: " ), expected.size() ), expected );
		}
#else
		}
#endif
	}
};  // class ThreeLetterCodeAmbiguityTests
