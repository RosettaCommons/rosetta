// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/forge/build/SegmentSwap.cxxtest.hh
/// @brief  unit tests for SegmentSwap BuildInstruction
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <protocols/forge/build/SegmentSwap.hh>

#include <numeric/xyzVector.hh>

#include <string>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


class SegmentSwapTests : public CxxTest::TestSuite
{


public: // setup


	typedef std::string String;
	typedef core::Size Size;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef core::kinematics::Edge Edge;
	typedef core::kinematics::FoldTree FoldTree;
	typedef core::kinematics::MoveMap MoveMap;
	typedef protocols::forge::build::Interval Interval;
	typedef protocols::forge::build::SegmentSwap SegmentSwap;


	SegmentSwapTests() {};


	// Shared initialization.
	void setUp() {
		core_init();
	}


	// Shared finalization.
	void tearDown() {
	}


public: // re-used methods


	/// @brief return a Pose with a continuous topology
	Pose continuous20_pose() {
		Pose pose;
		core::pose::make_pose_from_sequence(
			pose,
			"G[GLY:NtermProteinFull]FFFFFFFFFFFFFFFFFFG[GLY:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
		);

		for ( core::Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
			pose.set_secstruct( i, 'L' );
		}

		return pose;
	}


	/// @brief return a Pose with a continuous topology
	Pose continuous10_pose() {
		Pose pose;
		core::pose::make_pose_from_sequence(
			pose,
			"K[LYS:NtermProteinFull]EEEEEEEEK[LYS:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
		);

		for ( core::Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
			pose.set_secstruct( i, 'E' );
		}

		return pose;
	}


	/// @brief return a Pose with a cutpoint at 9 and jump from 7 to 14
	Pose cut_pose() {
		Pose pose;
		core::pose::make_pose_from_sequence(
			pose,
			"A[ALA:NtermProteinFull]CDEFGHVVVVVPQRSTVWY[TYR:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
		);

		for ( Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
			pose.set_secstruct( i, 'H' );
		}

		FoldTree ft;
		ft.add_edge( Edge( 1, 7, Edge::PEPTIDE ) );
		ft.add_edge( Edge( 7, 9, Edge::PEPTIDE ) );
		ft.add_edge( Edge( 10, 14, Edge::PEPTIDE ) );
		ft.add_edge( Edge( 14, 20, Edge::PEPTIDE ) );
		ft.add_edge( Edge( 7, 14, 1 ) ); // jump
		ft.reorder( 1 );

		pose.fold_tree( ft );

		for ( Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
			pose.set_phi( i, 10.0 * i );
			pose.set_psi( i, 5.0 *i );
			pose.set_omega( i, 180.0 );
		}

		pose.residue( 1 ); // force refold

		return pose;
	}


	/// @brief return a Pose with two cutpoints at 5 and 14 with jumps from 2->7 and 11->17
	Pose cut2_pose() {
		Pose pose;
		core::pose::make_pose_from_sequence(
			pose,
			"A[ALA:NtermProteinFull]WWWFGHIKLMNWWWSTVWY[TYR:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
		);

		for ( Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
			pose.set_secstruct( i, 'L' );
		}

		FoldTree ft;
		ft.simple_tree( 20 );
		ft.new_jump( 2, 7, 5 );
		ft.new_jump( 11, 17, 14 );

		pose.fold_tree( ft );

		for ( Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
			pose.set_phi( i, 5.0 * i );
			pose.set_psi( i, 10.0 *i );
			pose.set_omega( i, 180.0 );
		}

		pose.residue( 1 ); // force refold

		return pose;
	}


public: // tests


	/// @brief test simple swap of continuous segments
	void test_continuous_swap() {
		// create dummy poses
		Pose c20 = continuous20_pose();
		Pose c10 = continuous10_pose();

		// Store original position of CA of first, middle, last residue of cut2.
		// These shouldn't change after the operation.
		Vector first_CA = c10.residue( 1 ).xyz( "CA" );
		Vector middle_CA = c10.residue( 5 ).xyz( "CA" );
		Vector last_CA = c10.residue( 10 ).xyz( "CA" );

		// MoveMap for placing jumps
		MoveMap movemap; // empty, place jump anywhere

		// replace section of c20 with c10
		SegmentSwap swap( Interval( 4, 11 ), movemap, c10 );
		swap.modify( c20 );

		// artificially invoke change
		c20.set_phi( 2, 93.0 );
		c20.set_phi( 21, 39.0 );
		c20.residue( 2 ); // force coordinate update

		TS_ASSERT_EQUALS( swap.interval().left, 4 );
		TS_ASSERT_EQUALS( swap.interval().right, 13 );
		TS_ASSERT( swap.original_interval_valid() );
		TS_ASSERT_EQUALS( c20.n_residue(), 22 );
		TS_ASSERT_EQUALS( c20.fold_tree().num_cutpoint(), 2 );
		TS_ASSERT_EQUALS( c20.annotated_sequence(), "G[GLY:NtermProteinFull]FFKEEEEEEEEKFFFFFFFFG[GLY:CtermProteinFull]" );
		TS_ASSERT_EQUALS( c20.secstruct(), "LLLEEEEEEEEEELLLLLLLLL" );
		TS_ASSERT_EQUALS( c20.residue( 4 ).xyz( "CA" ), first_CA );
		TS_ASSERT_EQUALS( c20.residue( 8 ).xyz( "CA" ), middle_CA );
		TS_ASSERT_EQUALS( c20.residue( 13 ).xyz( "CA" ), last_CA );
	}


	/// @brief test swap of complex topologies, continuous -> complex, known jumps
	void test_continuous_complex_swap1() {
		// create dummy poses
		Pose cut = cut_pose();
		Pose cut2 = cut2_pose();

		// Store original position of CA of first, middle, last residue of cut2.
		// These shouldn't change after the operation.
		Vector first_CA = cut2.residue( 1 ).xyz( "CA" );
		Vector middle_CA = cut2.residue( 10 ).xyz( "CA" );
		Vector last_CA = cut2.residue( 20 ).xyz( "CA" );

		// MoveMap for placing jumps
		MoveMap movemap;
		movemap.set_bb_true_range( 1, cut.n_residue() );
		movemap.set_bb( 19 , false );

		// Replace section of cut with cut2.  The interval from [13, 17]
		// will blow away the original jump.  If working properly, the
		// procedure should connect all disconnected sections.
		SegmentSwap swap( Interval( 13, 17 ), movemap, cut2 );
		swap.modify( cut );

		TS_ASSERT_EQUALS( swap.interval().left, 13 );
		TS_ASSERT_EQUALS( swap.interval().right, 32 );
		TS_ASSERT( swap.original_interval_valid() );
		TS_ASSERT_EQUALS( cut.n_residue(), 35 );
		TS_ASSERT_EQUALS( cut.fold_tree().num_cutpoint(), 5 );
		TS_ASSERT_EQUALS( cut.annotated_sequence(), "A[ALA:NtermProteinFull]CDEFGHVVVVVAWWWFGHIKLMNWWWSTVWYVWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( cut.secstruct(), "HHHHHHHHHHHHLLLLLLLLLLLLLLLLLLLLHHH" );
		TS_ASSERT_EQUALS( cut.residue( 13 ).xyz( "CA" ), first_CA );
		TS_ASSERT_EQUALS( cut.residue( 22 ).xyz( "CA" ), middle_CA );
		TS_ASSERT_EQUALS( cut.residue( 32 ).xyz( "CA" ), last_CA );
		TS_ASSERT( cut.fold_tree().is_jump_point( 34 ) );
	}


	/// @brief test swap of complex topologies, continuous -> complex, random jumps
	void test_continuous_complex_swap2() {
		// create dummy poses
		Pose cut = cut_pose();
		Pose cut2 = cut2_pose();

		// Store original position of CA of first, middle, last residue of cut2.
		// These shouldn't change after the operation.
		Vector first_CA = cut2.residue( 1 ).xyz( "CA" );
		Vector middle_CA = cut2.residue( 10 ).xyz( "CA" );
		Vector last_CA = cut2.residue( 20 ).xyz( "CA" );

		// MoveMap for placing jumps
		MoveMap movemap;
		movemap.set_bb_true_range( 1, cut.n_residue() );
		movemap.set_bb( 18 , false );
		movemap.set_bb( 19 , false );
		movemap.set_bb( 20 , false );

		// Replace section of cut with cut2.  The interval from [12, 17]
		// will blow away the original jump.  If working properly, the
		// procedure should connect all disconnected sections.
		SegmentSwap swap( Interval( 13, 17 ), movemap, cut2 );
		swap.modify( cut );

		TS_ASSERT_EQUALS( swap.interval().left, 13 );
		TS_ASSERT_EQUALS( swap.interval().right, 32 );
		TS_ASSERT( swap.original_interval_valid() );
		TS_ASSERT_EQUALS( cut.n_residue(), 35 );
		TS_ASSERT_EQUALS( cut.fold_tree().num_cutpoint(), 5 );
		TS_ASSERT_EQUALS( cut.annotated_sequence(), "A[ALA:NtermProteinFull]CDEFGHVVVVVAWWWFGHIKLMNWWWSTVWYVWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( cut.secstruct(), "HHHHHHHHHHHHLLLLLLLLLLLLLLLLLLLLHHH" );
		TS_ASSERT_EQUALS( cut.residue( 13 ).xyz( "CA" ), first_CA );
		TS_ASSERT_EQUALS( cut.residue( 22 ).xyz( "CA" ), middle_CA );
		TS_ASSERT_EQUALS( cut.residue( 32 ).xyz( "CA" ), last_CA );
		TS_ASSERT( cut.fold_tree().is_jump_point( 33 ) || cut.fold_tree().is_jump_point( 34 ) || cut.fold_tree().is_jump_point( 35 ) );
	}

};
