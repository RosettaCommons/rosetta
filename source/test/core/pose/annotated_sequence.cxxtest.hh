// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/core/pose/annotated_sequence.cxxtest.hh
/// @brief   Test suite for utility functions for creating poses from annotated sequences.
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/pose/annotated_sequence.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

// Utility header
#include <utility/vector1.hh>


class AnnotatedSequenceTests : public CxxTest::TestSuite {
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
	// Confirm that an IUPAC sequence is correctly read and ResidueTypes correctly selected.
	void test_residue_types_from_saccharide_sequence()
	{
		using namespace utility;
		using namespace core;
		using namespace chemical;
		using namespace pose;

		ResidueTypeSetCOP residue_set( ChemicalManager::get_instance()->residue_type_set( "fa_standard" ) );
		vector1< ResidueTypeCOP > residue_types;

		TS_TRACE( "Testing a maltotriose, a simple linear sugar, relying on default settings." );
		residue_types = residue_types_from_saccharide_sequence(
				"Glcp-Glcp-Glcp", *residue_set );

		TS_ASSERT_EQUALS( residue_types.size(), 3 );
		TS_ASSERT_EQUALS( residue_types[ 1 ]->name(), "->4)-alpha-D-Glcp:non-reducing_end" );
		TS_ASSERT_EQUALS( residue_types[ 2 ]->name(), "->4)-alpha-D-Glcp" );
		TS_ASSERT_EQUALS( residue_types[ 3 ]->name(), "->4)-alpha-D-Glcp:branch_lower_terminus" );


		TS_TRACE( "Testing Lewisx, which has one branch and modified sugars." );
		residue_types = residue_types_from_saccharide_sequence(
				"beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpNAc", *residue_set );

		TS_ASSERT_EQUALS( residue_types.size(), 3 );
		TS_ASSERT_EQUALS( residue_types[ 1 ]->name(), "->4)-beta-D-Galp:non-reducing_end" );
		TS_ASSERT_EQUALS( residue_types[ 2 ]->name(), "->4)-alpha-L-Fucp:non-reducing_end:branch_lower_terminus" );
		TS_ASSERT_EQUALS( residue_types[ 3 ]->name(), "->4)-alpha-D-Glcp:branch_lower_terminus:->3)-branch:2-AcNH" );


		TS_TRACE( "Testing a 14-mer, which has nested branches." );
		residue_types = residue_types_from_saccharide_sequence(
				"a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Manp-(1->2)-a-D-Manp-(1->2)-a-D-Manp-(1->3)-"
				"[a-D-Manp-(1->2)-a-D-Manp-(1->3)-[a-D-Manp-(1->2)-a-D-Manp-(1->6)]-a-D-Manp-(1->6)]-b-D-Manp-(1->4)-"
				"b-D-GlcpNAc-(1->4)-b-D-GlcpNAc-", *residue_set );

		TS_ASSERT_EQUALS( residue_types.size(), 14 );
		TS_ASSERT_EQUALS( residue_types[ 1 ]->name(), "->4)-alpha-D-Glcp:non-reducing_end" );
		TS_ASSERT_EQUALS( residue_types[ 7 ]->name(), "->4)-alpha-D-Manp:non-reducing_end" );
		TS_ASSERT_EQUALS( residue_types[ 8 ]->name(), "->2)-alpha-D-Manp" );
		TS_ASSERT_EQUALS( residue_types[ 9 ]->name(), "->4)-alpha-D-Manp:non-reducing_end" );
		TS_ASSERT_EQUALS( residue_types[ 10 ]->name(), "->2)-alpha-D-Manp:branch_lower_terminus" );
		TS_ASSERT_EQUALS( residue_types[ 11 ]->name(), "->3)-alpha-D-Manp:branch_lower_terminus:->6)-branch" );
		TS_ASSERT_EQUALS( residue_types[ 12 ]->name(), "->3)-beta-D-Manp:->6)-branch" );
		TS_ASSERT_EQUALS( residue_types[ 14 ]->name(), "->4)-beta-D-Glcp:branch_lower_terminus:2-AcNH" );
	}

	// TODO: Somebody needs to add tests for the other functions in annotated_sequence.hh.
};  // class AnnotatedSequenceTests
