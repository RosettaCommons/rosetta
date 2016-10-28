// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/protocols/enzymatic_movers/GlycosyltransferaseMover.cxxtest.hh
/// @brief   Test suite for basic GlycosyltransferaseMover functionality.
/// @author  Labonte <JWLabonte@jhu.edu>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <protocols/enzymatic_movers/GlycosyltransferaseMover.hh>

// Project headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>

// Utility header
#include <utility/vector1.hh>


class GlycosyltransferaseMoverTests : public CxxTest::TestSuite {
public: // Standard methods ///////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init_with_additional_options( "-include_sugars" );
	}

	// Destruction
	void tearDown()
	{}


public: // Tests //////////////////////////////////////////////////////////////
	// Confirm that a simple pose can be glycosylated with the Mover.
	void test_apply()
	{
		using namespace core::pose;

		TS_TRACE( "Testing GlycosyltransferaseMover::apply()." );

		Pose peptide;
		make_pose_from_sequence( peptide, "NASANASA", "fa_standard" );

		TS_ASSERT_EQUALS( peptide.sequence(), "NASANASA" );
		TS_ASSERT_EQUALS( peptide.fold_tree().to_string(), "FOLD_TREE  EDGE 1 8 -1 " );

		protocols::enzymatic_movers::GlycosyltransferaseMover generic_N_glycosyltransferase;
		generic_N_glycosyltransferase.exclude_site( 5 );  // Limit to single site for testing purposes.
		generic_N_glycosyltransferase.ensure_site( 1 );  // Ensure it glycosylates for testing purposes.
		generic_N_glycosyltransferase.perform_major_reaction_only();  // can't have mixed results in a unit test!

		generic_N_glycosyltransferase.apply( peptide );

		TS_ASSERT_EQUALS( peptide.sequence(), "NASANASAZZZZZ" );
		TS_ASSERT_EQUALS( peptide.fold_tree().to_string(),
				"FOLD_TREE  EDGE 1 8 -1  "
				"EDGE 1 9 -2  ND2  C1   EDGE 9 12 -1  EDGE 11 13 -2  O6   C1  " );

		utility::vector1< core::uint > const excluded_sites( 1, 1 );
		utility::vector1< core::uint > const ensured_sites( 1, 5 );
		generic_N_glycosyltransferase.set_excluded_sites( excluded_sites );  // Change to a different site.
		generic_N_glycosyltransferase.set_ensured_sites( ensured_sites );

		generic_N_glycosyltransferase.apply( peptide );

		TS_ASSERT_EQUALS( peptide.sequence(), "NASANASAZZZZZZZZZZ" );
		TS_ASSERT_EQUALS( peptide.fold_tree().to_string(),
				"FOLD_TREE  EDGE 1 8 -1  "
				"EDGE 1 9 -2  ND2  C1   EDGE 9 12 -1  EDGE 11 13 -2  O6   C1   "
				"EDGE 5 14 -2  ND2  C1   EDGE 14 17 -1  EDGE 16 18 -2  O6   C1  ");
	}
};  // class GlycosyltransferaseMoverTests
