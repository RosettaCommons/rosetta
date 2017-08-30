// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>

// Utility header
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("core.pose.annotated_sequence.cxxtest");

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

		TR <<  "Testing a maltotriose, a simple linear sugar, relying on default settings."  << std::endl;
		residue_types = residue_types_from_saccharide_sequence(
			"Glcp-Glcp-Glcp", *residue_set );
		// Order 3-2-1

		TS_ASSERT_EQUALS( residue_types.size(), 3 );
		TS_ASSERT_EQUALS( residue_types[ 3 ]->name(), "->4)-alpha-D-Glcp:non-reducing_end" );
		TS_ASSERT_EQUALS( residue_types[ 2 ]->name(), "->4)-alpha-D-Glcp" );
		TS_ASSERT_EQUALS( residue_types[ 1 ]->name(), "->4)-alpha-D-Glcp" );


		TR <<  "Testing Lewisx, which has one branch and modified sugars."  << std::endl;
		residue_types = residue_types_from_saccharide_sequence(
			"beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpNAc", *residue_set );
		// Order 2-[3]-1

		TS_ASSERT_EQUALS( residue_types.size(), 3 );
		TS_ASSERT_EQUALS( residue_types[ 2 ]->name(), "->4)-beta-D-Galp:non-reducing_end" );
		TS_ASSERT_EQUALS( residue_types[ 3 ]->name(), "->4)-alpha-L-Fucp:non-reducing_end" );
		TS_ASSERT_EQUALS( residue_types[ 1 ]->name(), "->4)-alpha-D-Glcp:->3)-branch:2-AcNH" );


		TR <<  "Testing a 14-mer, which has nested branches."  << std::endl;
		residue_types = residue_types_from_saccharide_sequence(
			"a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Manp-(1->2)-a-D-Manp-(1->2)-a-D-Manp-(1->3)-"
			"[a-D-Manp-(1->2)-a-D-Manp-(1->3)-[a-D-Manp-(1->2)-a-D-Manp-(1->6)]-a-D-Manp-(1->6)]-b-D-Manp-(1->4)-"
			"b-D-GlcpNAc-(1->4)-b-D-GlcpNAc-", *residue_set );
		// Order 9-8-7-6-5-4-[12-11-[14-13]-10]-3-2-1

		TS_ASSERT_EQUALS( residue_types.size(), 14 );
		TS_ASSERT_EQUALS( residue_types[ 9 ]->name(), "->4)-alpha-D-Glcp:non-reducing_end" );
		TS_ASSERT_EQUALS( residue_types[ 12 ]->name(), "->4)-alpha-D-Manp:non-reducing_end" );
		TS_ASSERT_EQUALS( residue_types[ 11 ]->name(), "->2)-alpha-D-Manp" );
		TS_ASSERT_EQUALS( residue_types[ 14 ]->name(), "->4)-alpha-D-Manp:non-reducing_end" );
		TS_ASSERT_EQUALS( residue_types[ 13 ]->name(), "->2)-alpha-D-Manp" );
		TS_ASSERT_EQUALS( residue_types[ 10 ]->name(), "->3)-alpha-D-Manp:->6)-branch" );
		TS_ASSERT_EQUALS( residue_types[ 3 ]->name(), "->3)-beta-D-Manp:->6)-branch" );
		TS_ASSERT_EQUALS( residue_types[ 1 ]->name(), "->4)-beta-D-Glcp:2-AcNH" );


		TR <<  "Testing multiple non-nested branches."  << std::endl;
		residue_types = residue_types_from_saccharide_sequence(
			"beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)-alpha-L-Fucp-(1->3)]-beta-D-Galp-(1->4)-[a-D-Manp-(1->2)-a-D-Manp-(1->2)]-beta-D-Galp-(1->4)",
			*residue_set );
		// Order 3-[7-6]-2-[5-4]-1
		TS_ASSERT_EQUALS( residue_types[ 1 ]->name(), "->4)-beta-D-Galp:->2)-branch" );
		TS_ASSERT_EQUALS( residue_types[ 2 ]->name(), "->4)-beta-D-Galp:->3)-branch" );
		TS_ASSERT_EQUALS( residue_types[ 3 ]->name(), "->4)-beta-D-Galp:non-reducing_end" );
		TS_ASSERT_EQUALS( residue_types[ 4 ]->name(), "->2)-alpha-D-Manp" );
		TS_ASSERT_EQUALS( residue_types[ 5 ]->name(), "->4)-alpha-D-Manp:non-reducing_end" );
		TS_ASSERT_EQUALS( residue_types[ 6 ]->name(), "->3)-alpha-L-Fucp" );
		TS_ASSERT_EQUALS( residue_types[ 7 ]->name(), "->4)-alpha-L-Fucp:non-reducing_end" );
	}

	// Test that pose connectivity is set up correctly.
	void test_make_pose_from_saccharide_sequence()
	{
		core::pose::Pose pose;
		TR << "Testing making pose from simple sequence" << std::endl;
		make_pose_from_saccharide_sequence( pose, "Glcp-Glcp-Glcp" );

		TS_ASSERT_EQUALS( pose.residue(1).connected_residue_at_upper(), 2 );
		TS_ASSERT_EQUALS( pose.residue(2).connected_residue_at_upper(), 3 );
		TS_ASSERT_EQUALS( pose.residue(2).connected_residue_at_lower(), 1 );
		TS_ASSERT_EQUALS( pose.residue(3).connected_residue_at_lower(), 2 );

		TR <<  "Testing making pose from Lewisx, which has one branch and modified sugars."  << std::endl;
		make_pose_from_saccharide_sequence( pose, "beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpNAc" );
		// Should be constructed as 2-[3]-1
		TS_ASSERT_EQUALS( pose.residue_type(3).name(), "->4)-alpha-L-Fucp:non-reducing_end" );
		TS_ASSERT_EQUALS( pose.residue_type(2).name(), "->4)-beta-D-Galp:non-reducing_end" );
		TS_ASSERT_EQUALS( pose.residue(1).connected_residue_at_upper(), 2 );
		TS_ASSERT_EQUALS( pose.residue(2).connected_residue_at_lower(), 1 );
		TS_ASSERT_EQUALS( pose.residue(3).connected_residue_at_lower(), 1 );

		TR <<  "Testing making pose from a 14-mer, which has nested branches."  << std::endl;
		make_pose_from_saccharide_sequence( pose,
			"a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Manp-(1->2)-a-D-Manp-(1->2)-a-D-Manp-(1->3)-"
			"[a-D-Manp-(1->2)-a-D-Manp-(1->3)-[a-D-Manp-(1->2)-a-D-Manp-(1->6)]-a-D-Manp-(1->6)]"
			"-b-D-Manp-(1->4)-b-D-GlcpNAc-(1->4)-b-D-GlcpNAc-" );

		// Should be constructed as  9-8-7-6-5-4-[12-11-[14-13]-10]-3-2-1
		TS_ASSERT_EQUALS( pose.size(), 14 );
		TS_ASSERT_EQUALS( pose.residue_type(3).name(), "->3)-beta-D-Manp:->6)-branch" );
		TS_ASSERT_EQUALS( pose.residue_type(10).name(), "->3)-alpha-D-Manp:->6)-branch" );
		TS_ASSERT_EQUALS( pose.residue_type(14).name(), "->4)-alpha-D-Manp:non-reducing_end" );
		TS_ASSERT_EQUALS( pose.residue_type(12).name(), "->4)-alpha-D-Manp:non-reducing_end" );
		TS_ASSERT_EQUALS( pose.residue_type(9).name(), "->4)-alpha-D-Glcp:non-reducing_end" );

		TS_ASSERT_EQUALS( pose.residue(2).connected_residue_at_lower(), 1 );
		TS_ASSERT_EQUALS( pose.residue(3).connected_residue_at_lower(), 2 );
		TS_ASSERT_EQUALS( pose.residue(4).connected_residue_at_lower(), 3 );
		TS_ASSERT_EQUALS( pose.residue(5).connected_residue_at_lower(), 4 );
		TS_ASSERT_EQUALS( pose.residue(6).connected_residue_at_lower(), 5 );
		TS_ASSERT_EQUALS( pose.residue(7).connected_residue_at_lower(), 6 );
		TS_ASSERT_EQUALS( pose.residue(8).connected_residue_at_lower(), 7 );
		TS_ASSERT_EQUALS( pose.residue(9).connected_residue_at_lower(), 8 );
		TS_ASSERT_EQUALS( pose.residue(10).connected_residue_at_lower(), 3 );
		TS_ASSERT_EQUALS( pose.residue(11).connected_residue_at_lower(), 10 );
		TS_ASSERT_EQUALS( pose.residue(12).connected_residue_at_lower(), 11 );
		TS_ASSERT_EQUALS( pose.residue(13).connected_residue_at_lower(), 10 );
		TS_ASSERT_EQUALS( pose.residue(14).connected_residue_at_lower(), 13 );

		TS_ASSERT_EQUALS( pose.residue(1).connected_residue_at_upper(), 2 );
		TS_ASSERT_EQUALS( pose.residue(2).connected_residue_at_upper(), 3 );
		TS_ASSERT_EQUALS( pose.residue(3).connected_residue_at_upper(), 4 );
		TS_ASSERT_EQUALS( pose.residue(4).connected_residue_at_upper(), 5 );
		TS_ASSERT_EQUALS( pose.residue(5).connected_residue_at_upper(), 6 );
		TS_ASSERT_EQUALS( pose.residue(6).connected_residue_at_upper(), 7 );
		TS_ASSERT_EQUALS( pose.residue(7).connected_residue_at_upper(), 8 );
		TS_ASSERT_EQUALS( pose.residue(8).connected_residue_at_upper(), 9 );
		//TS_ASSERT_EQUALS( pose.residue(9).connected_residue_at_upper(), -- ); //
		TS_ASSERT_EQUALS( pose.residue(10).connected_residue_at_upper(), 11 );
		TS_ASSERT_EQUALS( pose.residue(11).connected_residue_at_upper(), 12 );
		//TS_ASSERT_EQUALS( pose.residue(12).connected_residue_at_upper(), -- ); //
		TS_ASSERT_EQUALS( pose.residue(13).connected_residue_at_upper(), 14 );

	}

	// TODO: Somebody needs to add tests for the other functions in annotated_sequence.hh.
};  // class AnnotatedSequenceTests
