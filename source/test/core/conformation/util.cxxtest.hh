// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/chemical/rings/util.cxxtest.hh
/// @brief   Test suite for conformation-related utility functions
/// @author  Labonte <JWLabonte@jhu.edu>


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Header
#include <core/conformation/util.hh>

// Package Header
#include <core/conformation/Residue.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

// Utility Header
#include <utility/vector1.hh>

// Basic Header
#include <basic/Tracer.hh>


static THREAD_LOCAL basic::Tracer TR( "core.conformation.util.cxxtest.hh" );


class ConformationUtilityFunctionTests : public CxxTest::TestSuite {
public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init_with_additional_options( "-include_sugars" );
	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////
	// Confirm that the position of a substituent on a ring can be correctly assigned.
	void test_position_of_atom_on_ring()
	{
		using namespace core::conformation;
		using namespace core::pose;

		TR << "Testing that that position_of_atom_on_ring() correctly figures out attachment positions." << std::endl;

		Pose monosaccharide;
		make_pose_from_saccharide_sequence( monosaccharide, "alpha-D-Glcp", "fa_standard" );
		Residue const & glucose( monosaccharide.residue( 1 ) );
		utility::vector1< core::uint > const & ring_atoms( glucose.type().ring_atoms( 1 ) );

		uint const O1_index( glucose.atom_index( " O1 " ) );
		uint const O2_index( glucose.atom_index( " O2 " ) );
		uint const O3_index( glucose.atom_index( " O3 " ) );
		uint const O4_index( glucose.atom_index( " O4 " ) );
		uint const O5_index( glucose.atom_index( " O5 " ) );
		uint const O6_index( glucose.atom_index( " O6 " ) );

		TS_ASSERT_EQUALS( position_of_atom_on_ring( glucose, O1_index, ring_atoms ), 1 );
		TS_ASSERT_EQUALS( position_of_atom_on_ring( glucose, O2_index, ring_atoms ), 2 );
		TS_ASSERT_EQUALS( position_of_atom_on_ring( glucose, O3_index, ring_atoms ), 3 );
		TS_ASSERT_EQUALS( position_of_atom_on_ring( glucose, O4_index, ring_atoms ), 4 );
		TS_ASSERT_EQUALS( position_of_atom_on_ring( glucose, O5_index, ring_atoms ), 5 );  // O5 is in the ring itself.
		TS_ASSERT_EQUALS( position_of_atom_on_ring( glucose, O6_index, ring_atoms ), 0 );  // O6 is exocyclic.
	}
};  // class ConformationUtilityFunctionTests
