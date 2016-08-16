// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/forge/build/SegmentRebuild.cxxtest.hh
/// @brief  unit tests for SegmentRebuild BuildInstruction
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
#include <protocols/forge/build/SegmentRebuild.hh>

#include <numeric/xyzVector.hh>

#include <string>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


class SegmentRebuildTests : public CxxTest::TestSuite
{


public: // setup


	typedef std::string String;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef core::kinematics::Edge Edge;
	typedef core::kinematics::FoldTree FoldTree;
	typedef protocols::forge::build::Interval Interval;
	typedef protocols::forge::build::SegmentRebuild SegmentRebuild;


	SegmentRebuildTests() {};


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


	/// @brief return a Pose with a cutpoint at 9 and jump from 7 to 14
	Pose cut_pose() {
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
		ft.add_edge( Edge( 1, 7, Edge::PEPTIDE ) );
		ft.add_edge( Edge( 7, 9, Edge::PEPTIDE ) );
		ft.add_edge( Edge( 10, 14, Edge::PEPTIDE ) );
		ft.add_edge( Edge( 14, 20, Edge::PEPTIDE ) );
		ft.add_edge( Edge( 7, 14, 1 ) ); // jump
		ft.reorder( 1 );

		pose.fold_tree( ft );

		return pose;
	}


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


	/// @brief test rebuild into continuous segment, same length
	void test_continuous_segment_internal_identity() {
		// create dummy pose
		Pose pose = continuous_pose();

		// store original position of CA of takeoff and landing
		// these shouldn't change after the operation
		Vector takeoff_CA = pose.residue( 4 ).xyz( "CA" );
		Vector landing_CA = pose.residue( 19 ).xyz( "CA" );

		// go from 14-mer to 14-mer
		SegmentRebuild rebuild( Interval( 5, 18 ), String( 14, 'H' ), "CAAAAAAAAAAAAT" ); // 14-mer, full-atom
		rebuild.modify( pose );

		// artificially invoke change
		pose.set_phi( 5, 93.0 );
		pose.set_psi( 18, 93.0 );
		pose.residue( 5 ); // force coordinate update

		TS_ASSERT_EQUALS( rebuild.interval().left, 5 );
		TS_ASSERT_EQUALS( rebuild.interval().right, 18 );
		TS_ASSERT_EQUALS( rebuild.original_interval().left, 5 );
		TS_ASSERT_EQUALS( rebuild.original_interval().right, 18 );
		TS_ASSERT( rebuild.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 20 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CDECAAAAAAAAAAAATWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLLLHHHHHHHHHHHHHHLL" );
		TS_ASSERT_EQUALS( pose.residue( 4 ).xyz( "CA" ), takeoff_CA );
		TS_ASSERT_EQUALS( pose.residue( 19 ).xyz( "CA" ), landing_CA );
	}


	/// @brief test rebuild into continuous segment, contracting length
	void test_continuous_segment_internal_contract() {
		// create dummy pose
		Pose pose = continuous_pose();

		// store original position of CA of takeoff and landing
		// these shouldn't change after the operation
		Vector takeoff_CA = pose.residue( 4 ).xyz( "CA" );
		Vector landing_CA = pose.residue( 19 ).xyz( "CA" );

		// go from 14-mer to 7-mer
		SegmentRebuild rebuild( Interval( 5, 18 ), String( 7, 'H' ), "CAAAAAT" ); // 7-mer, full-atom
		rebuild.modify( pose );

		// artificially invoke change
		pose.set_phi( 5, 93.0 );
		pose.set_psi( 11, 93.0 );
		pose.residue( 5 ); // force coordinate update

		TS_ASSERT_EQUALS( rebuild.interval().left, 5 );
		TS_ASSERT_EQUALS( rebuild.interval().right, 11 );
		TS_ASSERT_EQUALS( rebuild.original_interval().left, 5 );
		TS_ASSERT_EQUALS( rebuild.original_interval().right, 18 );
		TS_ASSERT( rebuild.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 13 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CDECAAAAATWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLLLHHHHHHHLL" );
		TS_ASSERT_EQUALS( pose.residue( 4 ).xyz( "CA" ), takeoff_CA );
		TS_ASSERT_EQUALS( pose.residue( 12 ).xyz( "CA" ), landing_CA );
	}


	/// @brief test rebuild into continuous segment, expanding length
	void test_continuous_segment_internal_expand() {
		// create dummy pose
		Pose pose = continuous_pose();

		// store original position of CA of takeoff and landing
		// these shouldn't change after the operation
		Vector takeoff_CA = pose.residue( 5 ).xyz( "CA" );
		Vector landing_CA = pose.residue( 8 ).xyz( "CA" );

		// go from 2-mer to 9-mer
		SegmentRebuild rebuild( Interval( 6, 7 ), String( 9, 'H' ), "CAAAAAAAT" ); // 9-mer, full-atom
		rebuild.modify( pose );

		// artificially invoke change
		pose.set_phi( 6, 93.0 );
		pose.set_psi( 14, 93.0 );
		pose.residue( 6 ); // force coordinate update

		TS_ASSERT_EQUALS( rebuild.interval().left, 6 );
		TS_ASSERT_EQUALS( rebuild.interval().right, 14 );
		TS_ASSERT_EQUALS( rebuild.original_interval().left, 6 );
		TS_ASSERT_EQUALS( rebuild.original_interval().right, 7 );
		TS_ASSERT( rebuild.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 27 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CDEFCAAAAAAATIKLMNPQRSTVWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLLLLHHHHHHHHHLLLLLLLLLLLLL" );
		TS_ASSERT_EQUALS( pose.residue( 5 ).xyz( "CA" ), takeoff_CA );
		TS_ASSERT_EQUALS( pose.residue( 15 ).xyz( "CA" ), landing_CA );
	}


	/// @brief test rebuild into cut segment, same length
	void test_cut_segment_identity() {
		// create dummy pose
		Pose pose = cut_pose();

		// store original position of CA of takeoff and landing
		// these shouldn't change after the operation
		Vector takeoff_CA = pose.residue( 4 ).xyz( "CA" );
		Vector landing_CA = pose.residue( 19 ).xyz( "CA" );

		// go from 14-mer to 14-mer
		SegmentRebuild rebuild( Interval( 5, 18 ), String( 14, 'H' ), "CAAAAAAAAAAAAT" ); // 14-mer, full-atom
		rebuild.modify( pose );

		// artificially invoke change
		pose.set_phi( 5, 93.0 );
		pose.set_psi( 18, 93.0 );
		pose.residue( 5 ); // force coordinate update

		TS_ASSERT_EQUALS( rebuild.interval().left, 5 );
		TS_ASSERT_EQUALS( rebuild.interval().right, 18 );
		TS_ASSERT_EQUALS( rebuild.original_interval().left, 5 );
		TS_ASSERT_EQUALS( rebuild.original_interval().right, 18 );
		TS_ASSERT( rebuild.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 20 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CDECAAAAAAAAAAAATWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLLLHHHHHHHHHHHHHHLL" );
		TS_ASSERT_EQUALS( pose.residue( 4 ).xyz( "CA" ), takeoff_CA );
		TS_ASSERT_EQUALS( pose.residue( 19 ).xyz( "CA" ), landing_CA );
	}


	/// @brief test rebuild into cut segment, contracting length
	void test_cut_segment_contract() {
		// create dummy pose
		Pose pose = cut_pose();

		// store original position of CA of takeoff and landing
		// these shouldn't change after the operation
		Vector takeoff_CA = pose.residue( 4 ).xyz( "CA" );
		Vector landing_CA = pose.residue( 19 ).xyz( "CA" );

		// go from 14-mer to 7-mer
		SegmentRebuild rebuild( Interval( 5, 18 ), String( 7, 'H' ), "CAAAAAT" ); // 7-mer, full-atom
		rebuild.modify( pose );

		// artificially invoke change
		pose.set_phi( 5, 93.0 );
		pose.set_psi( 11, 93.0 );
		pose.residue( 5 ); // force coordinate update

		TS_ASSERT_EQUALS( rebuild.interval().left, 5 );
		TS_ASSERT_EQUALS( rebuild.interval().right, 11 );
		TS_ASSERT_EQUALS( rebuild.original_interval().left, 5 );
		TS_ASSERT_EQUALS( rebuild.original_interval().right, 18 );
		TS_ASSERT( rebuild.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 13 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CDECAAAAATWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLLLHHHHHHHLL" );
		TS_ASSERT_EQUALS( pose.residue( 4 ).xyz( "CA" ), takeoff_CA );
		TS_ASSERT_EQUALS( pose.residue( 12 ).xyz( "CA" ), landing_CA );
	}


	/// @brief test rebuild into continuous segment, expanding length
	void test_cut_segment_expand() {
		// create dummy pose
		Pose pose = cut_pose();

		// store original position of CA of takeoff and landing
		// these shouldn't change after the operation
		Vector takeoff_CA = pose.residue( 8 ).xyz( "CA" );
		Vector landing_CA = pose.residue( 11 ).xyz( "CA" );

		// go from 2-mer to 9-mer
		SegmentRebuild rebuild( Interval( 9, 10 ), String( 9, 'H' ), "CAAAAAAAT" ); // 9-mer, full-atom
		rebuild.modify( pose );

		// artificially invoke change
		pose.set_phi( 9, 93.0 );
		pose.set_psi( 17, 93.0 );
		pose.residue( 9 ); // force coordinate update

		TS_ASSERT_EQUALS( rebuild.interval().left, 9 );
		TS_ASSERT_EQUALS( rebuild.interval().right, 17 );
		TS_ASSERT_EQUALS( rebuild.original_interval().left, 9 );
		TS_ASSERT_EQUALS( rebuild.original_interval().right, 10 );
		TS_ASSERT( rebuild.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 27 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CDEFGHICAAAAAAATMNPQRSTVWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLLLLLLLHHHHHHHHHLLLLLLLLLL" );
		TS_ASSERT_EQUALS( pose.residue( 8 ).xyz( "CA" ), takeoff_CA );
		TS_ASSERT_EQUALS( pose.residue( 18 ).xyz( "CA" ), landing_CA );
	}


	/// @brief test rebuild into double cut segment, same length
	void test_cut2_segment_identity() {
		// create dummy pose
		Pose pose = cut2_pose();

		// store original position of CA of takeoff and landing
		// these shouldn't change after the operation
		Vector takeoff_CA = pose.residue( 2 ).xyz( "CA" );
		Vector landing_CA = pose.residue( 19 ).xyz( "CA" );

		// go from 16-mer to 16-mer
		SegmentRebuild rebuild( Interval( 3, 18 ), String( 16, 'H' ), "CAAAAAAAAAAAAAAT" ); // 16-mer, full-atom
		rebuild.modify( pose );

		// artificially invoke change
		pose.set_phi( 3, 93.0 );
		pose.set_psi( 18, 93.0 );
		pose.residue( 3 ); // force coordinate update

		TS_ASSERT_EQUALS( rebuild.interval().left, 3 );
		TS_ASSERT_EQUALS( rebuild.interval().right, 18 );
		TS_ASSERT_EQUALS( rebuild.original_interval().left, 3 );
		TS_ASSERT_EQUALS( rebuild.original_interval().right, 18 );
		TS_ASSERT( rebuild.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 20 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CCAAAAAAAAAAAAAATWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLHHHHHHHHHHHHHHHHLL" );
		TS_ASSERT_EQUALS( pose.residue( 2 ).xyz( "CA" ), takeoff_CA );
		TS_ASSERT_EQUALS( pose.residue( 19 ).xyz( "CA" ), landing_CA );
	}


	/// @brief test rebuild into double cut segment, contract length
	void test_cut2_segment_contract() {
		// create dummy pose
		Pose pose = cut2_pose();

		// store original position of CA of takeoff and landing
		// these shouldn't change after the operation
		Vector takeoff_CA = pose.residue( 2 ).xyz( "CA" );
		Vector landing_CA = pose.residue( 19 ).xyz( "CA" );

		// go from 16-mer to 10-mer
		SegmentRebuild rebuild( Interval( 3, 18 ), String( 10, 'H' ), "CAAAAAAAAT" ); // 10-mer, full-atom
		rebuild.modify( pose );

		// artificially invoke change
		pose.set_phi( 3, 93.0 );
		pose.set_psi( 12, 93.0 );
		pose.residue( 3 ); // force coordinate update

		TS_ASSERT_EQUALS( rebuild.interval().left, 3 );
		TS_ASSERT_EQUALS( rebuild.interval().right, 12 );
		TS_ASSERT_EQUALS( rebuild.original_interval().left, 3 );
		TS_ASSERT_EQUALS( rebuild.original_interval().right, 18 );
		TS_ASSERT( rebuild.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 14 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CCAAAAAAAATWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLHHHHHHHHHHLL" );
		TS_ASSERT_EQUALS( pose.residue( 2 ).xyz( "CA" ), takeoff_CA );
		TS_ASSERT_EQUALS( pose.residue( 13 ).xyz( "CA" ), landing_CA );
	}


	/// @brief test rebuild into double cut segment, expand length
	void test_cut2_segment_expand() {
		// create dummy pose
		Pose pose = cut2_pose();

		// store original position of CA of takeoff and landing
		// these shouldn't change after the operation
		Vector takeoff_CA = pose.residue( 2 ).xyz( "CA" );
		Vector landing_CA = pose.residue( 19 ).xyz( "CA" );

		// go from 16-mer to 20-mer
		SegmentRebuild rebuild( Interval( 3, 18 ), String( 20, 'H' ), "CAAAAAAAAAAAAAAAAAAT" ); // 20-mer, full-atom
		rebuild.modify( pose );

		// artificially invoke change
		pose.set_phi( 3, 93.0 );
		pose.set_psi( 22, 93.0 );
		pose.residue( 3 ); // force coordinate update

		TS_ASSERT_EQUALS( rebuild.interval().left, 3 );
		TS_ASSERT_EQUALS( rebuild.interval().right, 22 );
		TS_ASSERT_EQUALS( rebuild.original_interval().left, 3 );
		TS_ASSERT_EQUALS( rebuild.original_interval().right, 18 );
		TS_ASSERT( rebuild.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 24 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CCAAAAAAAAAAAAAAAAAATWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLHHHHHHHHHHHHHHHHHHHHLL" );
		TS_ASSERT_EQUALS( pose.residue( 2 ).xyz( "CA" ), takeoff_CA );
		TS_ASSERT_EQUALS( pose.residue( 23 ).xyz( "CA" ), landing_CA );
	}


	/// @brief test rebuild of n-term of first chain of a two-chain segment
	void test_multi_chain_nterm_chain1() {
		// create dummy pose
		Pose pose = two_chain_pose();

		// store original position of CA of anchor
		// this shouldn't change after the operation
		Vector anchor_CA = pose.residue( 4 ).xyz( "CA" );

		// go from 3-mer to 5-mer
		SegmentRebuild rebuild( Interval( 1, 3 ), String( 5, 'H' ), "CAAAT" ); // 5-mer, full-atom
		rebuild.modify( pose );

		// artificially invoke change
		pose.set_psi( 5, 93.0 );
		pose.residue( 5 ); // force coordinate update

		TS_ASSERT_EQUALS( rebuild.interval().left, 1 );
		TS_ASSERT_EQUALS( rebuild.interval().right, 5 );
		TS_ASSERT_EQUALS( rebuild.original_interval().left, 1 );
		TS_ASSERT_EQUALS( rebuild.original_interval().right, 3 );
		TS_ASSERT( rebuild.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 22 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "C[CYS:NtermProteinFull]AAATEFGHIK[LYS:CtermProteinFull]L[LEU:NtermProteinFull]MNPQRSTVWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "HHHHHLLLLLLLLLLLLLLLLL" );
		TS_ASSERT_EQUALS( pose.residue( 6 ).xyz( "CA" ), anchor_CA );
	}


	/// @brief test rebuild of c-term of first chain of a two-chain segment
	void test_multi_chain_cterm_chain1() {
		// create dummy pose
		Pose pose = two_chain_pose();

		// store original position of CA of anchor
		// this shouldn't change after the operation
		Vector anchor_CA = pose.residue( 6 ).xyz( "CA" );

		// go from 3-mer to 1-mer
		SegmentRebuild rebuild( Interval( 7, 9 ), String( 1, 'H' ), "A" ); // 1-mer, full-atom
		rebuild.modify( pose );

		// artificially invoke change
		pose.set_phi( 7, 93.0 );
		pose.residue( 7 ); // force coordinate update

		TS_ASSERT_EQUALS( rebuild.interval().left, 7 );
		TS_ASSERT_EQUALS( rebuild.interval().right, 7 );
		TS_ASSERT_EQUALS( rebuild.original_interval().left, 7 );
		TS_ASSERT_EQUALS( rebuild.original_interval().right, 9 );
		TS_ASSERT( rebuild.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 18 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CDEFGA[ALA:CtermProteinFull]L[LEU:NtermProteinFull]MNPQRSTVWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLLLLLHLLLLLLLLLLL" );
		TS_ASSERT_EQUALS( pose.residue( 6 ).xyz( "CA" ), anchor_CA );
	}


	/// @brief test rebuild of n-term of second chain of a two-chain segment
	void test_multi_chain_nterm_chain2() {
		// create dummy pose
		Pose pose = two_chain_pose();

		// store original position of CA of anchor
		// this shouldn't change after the operation
		Vector anchor_CA = pose.residue( 15 ).xyz( "CA" );

		// go from 5-mer to 2-mer
		SegmentRebuild rebuild( Interval( 10, 14), String( 2, 'H' ), "AA" ); // 2-mer, full-atom
		rebuild.modify( pose );

		// artificially invoke change
		pose.set_psi( 11, 93.0 );
		pose.residue( 5 ); // force coordinate update

		TS_ASSERT_EQUALS( rebuild.interval().left, 10 );
		TS_ASSERT_EQUALS( rebuild.interval().right, 11 );
		TS_ASSERT_EQUALS( rebuild.original_interval().left, 10 );
		TS_ASSERT_EQUALS( rebuild.original_interval().right, 14 );
		TS_ASSERT( rebuild.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 17 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CDEFGHIK[LYS:CtermProteinFull]A[ALA:NtermProteinFull]ARSTVWY[TYR:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLLLLLLLLHHLLLLLL" );
		TS_ASSERT_EQUALS( pose.residue( 12 ).xyz( "CA" ), anchor_CA );
	}


	/// @brief test rebuild of c-term of second chain of a two-chain segment
	void test_multi_chain_cterm_chain2() {
		// create dummy pose
		Pose pose = two_chain_pose();

		// store original position of CA of anchor
		// this shouldn't change after the operation
		Vector anchor_CA = pose.residue( 16 ).xyz( "CA" );

		// identity, go from 4-mer to 4-mer
		SegmentRebuild rebuild( Interval( 17, 20 ), String( 4, 'H' ), "GAAG" ); // 4-mer, full-atom
		rebuild.modify( pose );

		// artificially invoke change
		pose.set_phi( 17, 93.0 );
		pose.residue( 17 ); // force coordinate update

		TS_ASSERT_EQUALS( rebuild.interval().left, 17 );
		TS_ASSERT_EQUALS( rebuild.interval().right, 20 );
		TS_ASSERT_EQUALS( rebuild.original_interval().left, 17 );
		TS_ASSERT_EQUALS( rebuild.original_interval().right, 20 );
		TS_ASSERT( rebuild.original_interval_valid() );
		TS_ASSERT_EQUALS( pose.n_residue(), 20 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( pose.annotated_sequence(), "A[ALA:NtermProteinFull]CDEFGHIK[LYS:CtermProteinFull]L[LEU:NtermProteinFull]MNPQRSGAAG[GLY:CtermProteinFull]" );
		TS_ASSERT_EQUALS( pose.secstruct(), "LLLLLLLLLLLLLLLLHHHH" );
		TS_ASSERT_EQUALS( pose.residue( 16 ).xyz( "CA" ), anchor_CA );
	}


};
