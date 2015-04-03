// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/forge/build/RelativeConnectRight.cxxtest.hh
/// @brief  unit tests for RelativeConnectRight BuildInstruction
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/RelativeConnectRight.hh>
#include <protocols/forge/build/RelativeSequencePosition.hh>
#include <protocols/forge/build/SegmentInsert.hh>

#include <string>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


class RelativeConnectRightTests : public CxxTest::TestSuite
{


public: // setup


	typedef std::string String;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef core::kinematics::Edge Edge;
	typedef core::kinematics::FoldTree FoldTree;
	typedef core::kinematics::Jump Jump;
	typedef protocols::forge::build::BuildManager BuildManager;
	typedef protocols::forge::build::CountFromLeft CountFromLeft;
	typedef protocols::forge::build::CountFromLeftOP CountFromLeftOP;
	typedef protocols::forge::build::Interval Interval;
	typedef protocols::forge::build::RelativeConnectRight RelativeConnectRight;
	typedef protocols::forge::build::RelativeConnectRightOP RelativeConnectRightOP;
	typedef protocols::forge::build::RelativeSequencePosition RelativeSequencePosition;
	typedef protocols::forge::build::RelativeSequencePositionOP RelativeSequencePositionOP;
	typedef protocols::forge::build::SegmentInsert SegmentInsert;
	typedef protocols::forge::build::SegmentInsertOP SegmentInsertOP;


	RelativeConnectRightTests() {};


	// Shared initialization.
	void setUp() {
		core_init();
	}


	// Shared finalization.
	void tearDown() {
	}


public: // re-used methods


	/// @brief return a Pose with a cutpoint at 9 and jump from 7 to 14
	Pose cut_pose() {
		Pose pose;
		core::pose::make_pose_from_sequence(
			pose,
			"A[ALA:NtermProteinFull]CDEFGHIKLMNPQRSTVWY[TYR:CtermProteinFull]",
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
			"A[ALA:NtermProteinFull]WWWFGHIKLMNPQRSTVWY[TYR:CtermProteinFull]",
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


	/// @brief test connection
	void test_connection() {
		using protocols::forge::build::SegmentInsertConnectionScheme::C;

		// create dummy pose
		Pose cut = cut_pose(); // jump from 7 to 14
		Pose cut2 = cut2_pose();

		// split the 'cut' pose into left and right halves
		Pose left_cut;
		for ( Size i = 1; i <= 9; ++i ) {
			left_cut.append_residue_by_bond( cut.residue( i ) );
		}

		Pose right_cut;
		for ( Size i = 10, ie = cut.n_residue(); i <= ie; ++i ) {
			right_cut.append_residue_by_bond( cut.residue( i ) );
		}

		// Store the distances between a CA on left_cut and a CA on right_cut.
		// This shouldn't change after the operation.
		Real const ca5_ca14 = cut.residue( 5 ).xyz( "CA" ).distance( cut.residue( 14 ).xyz( "CA" ) );
		Real const ca2_ca17 = cut.residue( 2 ).xyz( "CA" ).distance( cut.residue( 17 ).xyz( "CA" ) );
		Real const ca8_ca19 = cut.residue( 8 ).xyz( "CA" ).distance( cut.residue( 19 ).xyz( "CA" ) );
		Real const ca2_ca11 = cut.residue( 2 ).xyz( "CA" ).distance( cut.residue( 11 ).xyz( "CA" ) );

		// make an independent instruction
		SegmentInsertOP si( new SegmentInsert( Interval( 13, 15 ), "L^L", left_cut, true, C ) );

		// Make a RelativeConnectRight that depends on the SegmentInsert and
		// set it up so that it tries to mirror the jump in 'cut'.  We'll be
		// testing to see if everything is successful by testing the relative
		// distances between 'left_cut' and 'right_cut' and making sure they're
		// the same after the modifications.
		CountFromLeftOP cfl( new CountFromLeft() );
		cfl->left_skip = 1;
		cfl->p = 7;

		RelativeConnectRightOP rcr( new RelativeConnectRight( cfl, 5, right_cut ) );

		// Set up the rt by grabbing it from cut_pose.
		rcr->extract_rt( cut, 7, 14 );
		rcr->use_rt( true );

		// setup instructions and dependencies in manager
		BuildManager manager;
		manager.add( si );
		manager.add( rcr );
		manager.create_directed_dependency( si, rcr );

		// do the modify
		manager.modify( cut2 );

		// force refold
		cut2.residue( 1 );

		// test for distances
		Real const test_5_14 = cut2.residue( 18 ).xyz( "CA" ).distance( cut2.residue( 33 ).xyz( "CA" ) );
		Real const test_2_17 = cut2.residue( 15 ).xyz( "CA" ).distance( cut2.residue( 36 ).xyz( "CA" ) );
		Real const test_8_19 = cut2.residue( 21 ).xyz( "CA" ).distance( cut2.residue( 38 ).xyz( "CA" ) );
		Real const test_2_11 = cut2.residue( 15 ).xyz( "CA" ).distance( cut2.residue( 30 ).xyz( "CA" ) );

		TS_ASSERT_DELTA( ca5_ca14, test_5_14, 0.001 );
		TS_ASSERT_DELTA( ca2_ca17, test_2_17, 0.001 );
		TS_ASSERT_DELTA( ca8_ca19, test_8_19, 0.001 );
		TS_ASSERT_DELTA( ca2_ca11, test_2_11, 0.001 );

	}

};
