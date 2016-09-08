// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/forge/build/ConnectRight.cxxtest.hh
/// @brief  unit tests for ConnectRight BuildInstruction
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
#include <protocols/forge/build/ConnectRight.hh>

#include <string>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


class ConnectRightTests : public CxxTest::TestSuite
{


public: // setup


	typedef std::string String;
	typedef core::Size Size;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef core::kinematics::Edge Edge;
	typedef core::kinematics::FoldTree FoldTree;
	typedef protocols::forge::build::Interval Interval;
	typedef protocols::forge::build::ConnectRight ConnectRight;


	ConnectRightTests() {};


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

		for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
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

		for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
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

		for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
			pose.set_secstruct( i, 'L' );
		}

		FoldTree ft;
		ft.simple_tree( 20 );
		ft.new_jump( 2, 7, 5 );
		ft.new_jump( 11, 17, 14 );

		pose.fold_tree( ft );

		for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
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
		// create dummy pose
		Pose cut = cut_pose();
		Pose cut2 = cut2_pose();

		// Store original position of CA of first, middle, last residue of cut2.
		// These shouldn't change after the operation.
		Vector first_CA = cut2.residue( 1 ).xyz( "CA" );
		Vector middle_CA = cut2.residue( 10 ).xyz( "CA" );
		Vector last_CA = cut2.residue( 20 ).xyz( "CA" );

		// connect cut2 to the right of cut
		ConnectRight connect( 6, 10, cut2 );
		connect.modify( cut );

		// artificially invoke change
		cut.set_phi( 9, 93.0 );
		cut.residue( 9 ); // force coordinate update

		TS_ASSERT_EQUALS( connect.interval().left, 21 );
		TS_ASSERT_EQUALS( connect.interval().right, 40 );
		TS_ASSERT( !connect.original_interval_valid() );
		TS_ASSERT_EQUALS( cut.size(), 40 );
		TS_ASSERT_EQUALS( cut.fold_tree().num_cutpoint(), 4 );
		TS_ASSERT_EQUALS( cut.annotated_sequence(), "A[ALA:NtermProteinFull]CDEFGHIKLMNPQRSTVWY[TYR:CtermProteinFull]A[ALA:NtermProteinFull]WWWFGHIKLMNPQRSTVWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( cut.secstruct(), "HHHHHHHHHHHHHHHHHHHHLLLLLLLLLLLLLLLLLLLL" );
		TS_ASSERT_EQUALS( cut.residue( 21 ).xyz( "CA" ), first_CA );
		TS_ASSERT_EQUALS( cut.residue( 30 ).xyz( "CA" ), middle_CA );
		TS_ASSERT_EQUALS( cut.residue( 40 ).xyz( "CA" ), last_CA );
	}

};
