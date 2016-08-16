// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/BetaTurnDetectionFeatures.cxxtest.hh
/// @brief test suite for protocols/features/BetaTurnDetectionFeatures.hh/cc
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <protocols/features/BetaTurnDetection.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>


class BetaTurnDetectionFeaturesTest : public CxxTest::TestSuite {

public:
	void setUp() {
		protocols_init();
	}

	void test_ramachandran_hashes() {
		using namespace protocols::features;
		BetaTurnDetectionOP beta_turns( new BetaTurnDetection );

		// Test edge cases and more normal bounds for trans peptide planes for the case of phi <= 0
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0., 50., 180. ) == "A" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( -150., 0., 180. ) == "A" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0., -99.9, 180. ) == "A" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0., 50.1, 180. ) == "B" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0., -100., 180. ) == "B" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( -150., -150., 180. ) == "B" );

		// Test edge cases and more normal bounds for cis peptide planes for the case of phi <= 0
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0., 50., 0. ) == "a" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( -150., 0., 0. ) == "a" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0., -99.9, 0. ) == "a" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0., 50.1, 0. ) == "b" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0., -100., 0. ) == "b" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( -150., -150., 0. ) == "b" );

		// Test edge cases and more normal bounds for trans peptide planes for the case of phi >= 0
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0.1, 100., 180. ) == "L" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 150., 0., 180. ) == "L" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0.1, -49.9, 180. ) == "L" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0.1, 100.1, 180. ) == "E" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0.1, -50., 180. ) == "E" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 150., -150., 180. ) == "E" );

		// Test edge cases and more normal bounds for cis peptide planes for the case of phi >= 0
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0.1, 100., 0. ) == "l" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 150., 0., 0. ) == "l" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0.1, -49.9, 0. ) == "l" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0.1, 100.1, 0. ) == "e" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 0.1, -50., 0. ) == "e" );
		TS_ASSERT( beta_turns->determine_ramachandran_hash_for_residue_with_dihedrals( 150., -150., 0. ) == "e" );

	}

	void test_validate_ramachandran_hash() {
		using namespace protocols::features;
		BetaTurnDetectionOP beta_turns( new BetaTurnDetection );

		std::string test_string;

		// Test all trans, first residue is A
		test_string = "AA";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "AA" );

		test_string = "AB";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "AB" );

		test_string = "AL";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "AL" );

		test_string = "AE";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "AE" );

		// Test all trans, first residue is B
		test_string = "BA";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "BA" );

		test_string = "BB";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "BB" );

		test_string = "BL";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "BL" );

		test_string = "BE";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "BE" );

		// Test all trans, first residue is L
		test_string = "LA";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "LA" );

		test_string = "LB";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "LB" );

		test_string = "LL";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "LL" );

		test_string = "LE";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "LE" );

		// Test all trans, first residue is E
		test_string = "EA";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "EA" );

		test_string = "EB";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "EB" );

		test_string = "EL";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "EL" );

		test_string = "EE";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "EE" );

		// Test the well characterized beta-turns with a cis residue at position 3
		test_string = "Ba";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "Ba" );

		test_string = "Bb";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "Bb" );

		// Test poorly characterized cis residues
		test_string = "Bl";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "Xx" );

		test_string = "bb";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "xx" );

		test_string = "bL";
		beta_turns->validate_ramachandran_hash( test_string );
		TS_ASSERT( test_string == "xX" );

		// Test something that should cause an exception to be thrown
		try {
			test_string = "XX";
			beta_turns->validate_ramachandran_hash( test_string );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_error_message = "The Ramachandran hash 'XX' contains 'X,' which is not valid. Valid Ramachandran hashes are 'A', 'B', 'L' and 'E' for trans peptide bonds, and 'a', 'b', 'l' and 'e' for cis peptide bonds.";
			TS_ASSERT( expected_error_message == e.msg() );
		}

	}

};
