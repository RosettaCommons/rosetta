// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/forge/build/BuildManager.cxxtest.hh
/// @brief  unit tests for GrowLeft BuildInstruction
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/GrowLeft.hh>
#include <protocols/forge/build/GrowRight.hh>
#include <protocols/forge/build/SegmentRebuild.hh>

#include <string>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


class BuildManagerTests : public CxxTest::TestSuite
{


public: // setup


	typedef std::string String;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef protocols::forge::build::BuildManager BuildManager;
	typedef protocols::forge::build::BuildInstructionOP BuildInstructionOP;
	typedef protocols::forge::build::GrowLeft GrowLeft;
	typedef protocols::forge::build::GrowRight GrowRight;
	typedef protocols::forge::build::Interval Interval;
	typedef protocols::forge::build::SegmentRebuild SegmentRebuild;


	BuildManagerTests() {};


	// Shared initialization.
	void setUp() {
		core_init();
	}


	// Shared finalization.
	void tearDown() {
	}


public: // re-used methods


	/// @brief return a Pose with a continuous topology
	Pose continuous_pose() {
		Pose pose;
		core::pose::make_pose_from_sequence(
			pose,
			"A[ALA:NtermProteinFull]CDEFGHIKLMNPQRSTVWY[TYR:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
		);

		for ( core::Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
			pose.set_secstruct( i, 'L' );
		}

		return pose;
	}


public: // tests


	/// @brief test to make sure compatibility check in working order
	void test_compatibility_check() {
		BuildManager manager;
		manager.add( new GrowLeft( 1, String( 3, 'H' ) ) );
		manager.add( new SegmentRebuild( Interval( 1, 10 ), String( 3, 'H' ) ) );
		TS_ASSERT( !manager.compatibility_check() );

		manager.clear();
		manager.add( new GrowLeft( 1, String( 3, 'H' ) ) );
		manager.add( new SegmentRebuild( Interval( 2, 10 ), String( 3, 'H' ) ) );
		TS_ASSERT( manager.compatibility_check() );

		manager.clear();
		manager.add( new GrowRight( 20, String( 3, 'H' ) ) );
		manager.add( new SegmentRebuild( Interval( 10, 20 ), String( 3, 'H' ) ) );
		TS_ASSERT( !manager.compatibility_check() );

		manager.clear();
		manager.add( new GrowRight( 20, String( 3, 'H' ) ) );
		manager.add( new SegmentRebuild( Interval( 10, 19 ), String( 3, 'H' ) ) );
		TS_ASSERT( manager.compatibility_check() );

		manager.clear();
		manager.add( new GrowRight( 18, String( 3, 'H' ) ) );
		manager.add( new GrowLeft( 12, String( 3, 'H' ) ) );
		manager.add( new SegmentRebuild( Interval( 10, 20 ), String( 3, 'H' ) ) );
		TS_ASSERT( !manager.compatibility_check() );

		manager.clear();
		manager.add( new GrowLeft( 1, String( 3, 'H' ) ) );
		manager.add( new SegmentRebuild( Interval( 2, 19 ), String( 3, 'H' ) ) );
		manager.add( new GrowRight( 20, String( 3, 'H' ) ) );
		TS_ASSERT( manager.compatibility_check() );
	}


	/// @brief test modify operation
	void test_modify() {
		Pose pose = continuous_pose();

		BuildManager manager;

		manager.add( new SegmentRebuild( Interval( 5, 9 ), String( 3, 'H' ) ) );
		manager.add( new GrowLeft( 1, String( 2, 'H' ) ) );
		manager.add( new SegmentRebuild( Interval( 12, 16 ), String( 7, 'H' ) ) );
		manager.add( new GrowRight( pose.n_residue(), String( 4, 'H' ) ) );

		// test position mapping doesn't exist prior to modify()
		TS_ASSERT( manager.original2modified().empty() );
		TS_ASSERT( manager.intervals().empty() );
		TS_ASSERT( manager.positions().empty() );

		manager.modify( pose );
		TS_ASSERT_EQUALS( pose.n_residue(), 26 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 2 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]AACDEAAALMAAAAAAATVWYAAAA[ALA:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "HHLLLLHHHLLHHHHHHHLLLLHHHH" );

		// test position mapping exists post modify()
		TS_ASSERT( !manager.original2modified().empty() );
		TS_ASSERT( !manager.intervals().empty() );
		TS_ASSERT( !manager.positions().empty() );
	}


	/// @brief test dependency operations
	void test_dependencies() {
		BuildManager manager;

		// dummy instructions
		BuildInstructionOP sr = new SegmentRebuild( Interval( 1, 10 ), String ( 3, 'H' ) );
		BuildInstructionOP gl = new GrowLeft( 1, String( 2, 'H' ) );
		BuildInstructionOP gr = new GrowRight( 20, String( 7, 'E' ) );

		manager.add( sr );
		manager.add( gl );
		manager.add( gr );

		manager.create_directed_dependency( sr, gl );
		manager.create_directed_dependency( sr, gr );

		// dependencies added correctly?
		TS_ASSERT( !sr->has_dependencies() );
		TS_ASSERT_EQUALS( gl->n_dependencies(), 1 );
		TS_ASSERT_EQUALS( gr->n_dependencies(), 1 );
		TS_ASSERT_EQUALS( manager.n_dependencies(), 2 );

		// do a copy construct and see if dependencies are copied
		BuildManager manager2( manager );
		BuildManager::BIOPConstIterator ii = manager2.begin();

		TS_ASSERT_EQUALS( manager2.n_dependencies(), 2 );
		TS_ASSERT_EQUALS( (**ii).n_dependencies(), 0 );
		++ii;
		TS_ASSERT_EQUALS( (**ii).n_dependencies(), 1 );
		++ii;
		TS_ASSERT_EQUALS( (**ii).n_dependencies(), 1 );
	}


};
