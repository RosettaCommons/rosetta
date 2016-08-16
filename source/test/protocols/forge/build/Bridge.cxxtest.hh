// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/forge/build/Bridge.cxxtest.hh
/// @brief  unit tests for Bridge BuildInstruction
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <protocols/forge/build/Bridge.hh>

#include <string>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


class BridgeTests : public CxxTest::TestSuite
{


public: // setup


	typedef std::string String;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef core::kinematics::Edge Edge;
	typedef core::kinematics::FoldTree FoldTree;
	typedef protocols::forge::build::Bridge Bridge;
	typedef protocols::forge::build::Interval Interval;


	BridgeTests() {};


	// Shared initialization.
	void setUp() {
		core_init();
	}


	// Shared finalization.
	void tearDown() {
	}


public: // re-used methods


	/// @brief return a Pose with two cutpoints at 5 and 14 with jumps from 2->7 and 11->17
	Pose cut2_pose() {
		Pose pose;
		core::pose::make_pose_from_sequence(
			pose,
			"A[ALA:NtermProteinFull]CDEFGHIKLMNPQRSTVWY[TYR:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
		);

		for ( core::Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
			pose.set_secstruct( i, 'L' );
		}

		FoldTree ft;
		ft.simple_tree( 20 );
		ft.new_jump( 2, 7, 5 );
		ft.new_jump( 11, 17, 14 );

		pose.fold_tree( ft );

		return pose;
	}


	/// @brief return a two-chain Pose ( 9 res + 11 res )
	Pose two_chain_pose() {
		Pose pose;
		core::pose::make_pose_from_sequence(
			pose,
			"A[ALA:NtermProteinFull]CDEFGHIK[LYS:CtermProteinFull]L[LEU:NtermProteinFull]MNPQRSTVWY[TYR:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
		);

		for ( core::Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
			pose.set_secstruct( i, 'L' );
		}

		return pose;
	}


public: // tests


	/// @brief test bridging with a 0 residue bridge
	void test_bridge_size_zero() {
		Pose pose = cut2_pose();

		Bridge bridge( Interval( 5, 6 ), String() );
		bridge.modify( pose );

		// artificially invoke change
		pose.set_phi( 5, 93.0 );
		pose.residue( 5 ); // force coordinate update

		TS_ASSERT_EQUALS( bridge.interval().left, 5 );
		TS_ASSERT_EQUALS( bridge.interval().right, 6 );
		TS_ASSERT_EQUALS( bridge.original_interval().left, 5 );
		TS_ASSERT_EQUALS( bridge.original_interval().right, 6 );
		TS_ASSERT( bridge.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 20 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT( !pose.fold_tree().is_cutpoint( 5 ) );
		TS_ASSERT( pose.fold_tree().is_cutpoint( 14 ) );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CDEFGHIKLMNPQRSTVWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), String( 20, 'L' ) );
	}


	/// @brief test bridging with a 3 residue bridge
	void test_bridge_size_three() {
		Pose pose = cut2_pose();

		Bridge bridge( Interval( 14, 15 ), "HHH", "YYY" );
		bridge.modify( pose );

		// artificially invoke change
		pose.set_phi( 14, 93.0 );
		pose.residue( 14 ); // force coordinate update

		TS_ASSERT_EQUALS( bridge.interval().left, 14 );
		TS_ASSERT_EQUALS( bridge.interval().right, 18 );
		TS_ASSERT_EQUALS( bridge.original_interval().left, 14 );
		TS_ASSERT_EQUALS( bridge.original_interval().right, 15 );
		TS_ASSERT( bridge.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 23 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT( pose.fold_tree().is_cutpoint( 5 ) );
		TS_ASSERT( !pose.fold_tree().is_cutpoint( 14 ) );
		TS_ASSERT( !pose.fold_tree().is_cutpoint( 17 ) ); // if cutpoint was kept, it would be shifted up 3
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CDEFGHIKLMNPQYYYRSTVWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLLLLLLLLLLLLLHHHLLLLLL" );
	}


	/// @brief test bridging a two chain pose to make sure termini variants
	///  are removed
	void test_bridge_two_chain() {
		Pose pose = two_chain_pose();

		Bridge bridge( Interval( 9, 10 ), "HHHHH", "WWWWW" );
		bridge.modify( pose );

		// artificially invoke change
		pose.set_phi( 9, 93.0 );
		pose.residue( 9 ); // force coordinate update

		TS_ASSERT_EQUALS( bridge.interval().left, 9 );
		TS_ASSERT_EQUALS( bridge.interval().right, 15 );
		TS_ASSERT_EQUALS( bridge.original_interval().left, 9 );
		TS_ASSERT_EQUALS( bridge.original_interval().right, 10 );
		TS_ASSERT( bridge.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 25 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 0 );
		TS_ASSERT( !pose.fold_tree().is_cutpoint( 9 ) );
		TS_ASSERT( !pose.fold_tree().is_cutpoint( 14 ) ); // if cutpoint was kept, it would be shifted up 5
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CDEFGHIKWWWWWLMNPQRSTVWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLLLLLLLLHHHHHLLLLLLLLLLL" );
	}


};
