// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/forge/build/SegmentInsert.cxxtest.hh
/// @brief  unit tests for SegmentInsert BuildInstruction
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <protocols/forge/build/SegmentInsert.hh>

#include <numeric/angle.functions.hh>

#include <cmath>
#include <string>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <core/kinematics/AtomPointer.hh>
#include <core/kinematics/DomainMap.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>
#include <utility/io/all.fwd.hh>
#include <utility/io/icstream.fwd.hh>
#include <utility/io/ocstream.fwd.hh>
#include <numeric/BodyPosition.fwd.hh>
#include <numeric/Quaternion.fwd.hh>
#include <numeric/all.fwd.hh>


class SegmentInsertTests : public CxxTest::TestSuite
{


public: // setup


	typedef std::string String;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef core::kinematics::Edge Edge;
	typedef core::kinematics::FoldTree FoldTree;
	typedef core::kinematics::MoveMap MoveMap;
	typedef protocols::forge::build::Interval Interval;
	typedef protocols::forge::build::SegmentInsert SegmentInsert;


	SegmentInsertTests() {};


	// Shared initialization.
	void setUp() {
		core_init();
	}


	// Shared finalization.
	void tearDown() {
	}


public: // re-used methods


	/// @brief return a Pose with a continuous topology, helical geometry, and a distorted bond
	Pose helix10_distorted_pose() {
		using core::id::AtomID;

		Pose pose;
		core::pose::make_pose_from_sequence(
			pose,
			"K[LYS:NtermProteinFull]EEEEEEEEK[LYS:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
		);

		for ( core::Size i = 1, ie = pose.size(); i <= ie; ++i ) {
			pose.set_secstruct( i, 'H' );
		}

		for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
			pose.set_phi( i, -60.0 );
			pose.set_psi( i, -45.0 );
			pose.set_omega( i, 180.0 );
		}

		// elongate the N -> CA bond at residue 7
		AtomID const res7_n( pose.residue( 7 ).mainchain_atom( 1 ), 7 );
		AtomID const res7_ca( pose.residue( 7 ).mainchain_atom( 2 ), 7 );
		pose.conformation().set_bond_length( res7_n, res7_ca, 5.0 );

		return pose;
	}


	/// @brief return a Pose with a continuous topology and nonsensical geometry
	Pose continuous10_pose() {
		Pose pose;
		core::pose::make_pose_from_sequence(
			pose,
			"R[ARG:NtermProteinFull]DDDDDDDDR[ARG:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
		);

		for ( core::Size i = 1, ie = pose.size(); i <= ie; ++i ) {
			pose.set_secstruct( i, 'L' );
		}

		for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
			pose.set_phi( i, 11.0 * i );
			pose.set_psi( i, 7.0 * i );
			pose.set_omega( i, 180.0 );
		}

		return pose;
	}


	/// @brief return an extended Pose with a cutpoint at 9 and jump from 7 to 14
	///  and nonsensical geometry
	Pose cut_pose() {
		Pose pose;
		core::pose::make_pose_from_sequence(
			pose,
			"A[ALA:NtermProteinFull]CDEFGHVVVVVPQRSTVWY[TYR:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
		);

		for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
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

		for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
			pose.set_phi( i, 10.0 * i );
			pose.set_psi( i, 5.0 *i );
			pose.set_omega( i, 180.0 );
		}

		pose.residue( 1 ); // force refold

		return pose;
	}


	/// @brief return a Pose with two cutpoints at 5 and 14 with jumps from 2->7 and 11->17
	///  and nonsensical geometry
	Pose cut2_pose() {
		Pose pose;
		core::pose::make_pose_from_sequence(
			pose,
			"A[ALA:NtermProteinFull]WWWWWWWWWWWWWWWWW WY[TYR:CtermProteinFull]",
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


	/// @brief test insertion with N-side connection, simple replacement, completely
	///  internal insertion
	void test_N_side_connection_simple_internal() {
		using namespace protocols::forge::build::SegmentInsertConnectionScheme;

		using numeric::principal_angle_degrees;
		using numeric::nonnegative_principal_angle_degrees;

		// create dummy poses
		Pose h10 = helix10_distorted_pose();
		Pose c10 = continuous10_pose();

		// Store CA <-> CA distances between first, middle, and last residues
		// of h10 pairwise.  These shouldn't change after the operation.
		Real const f_m_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 5 ).xyz( "CA" ) );
		Real const m_l_dist = h10.residue( 5 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );
		Real const f_l_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );

		// Store CA positions of flanking left-1 and right+1 residues.
		// These shouldn't change after the operation.
		Vector const left_minus_one_CA = c10.residue( 3 ).xyz( "CA" );
		Vector const right_plus_one_CA = c10.residue( 9 ).xyz( "CA" );

		// Store bb torsions at endpoints; we are requesting that SegmentInsert
		// retain them.
		Real const left_endpoint_phi = c10.phi( 4 );
		Real const right_endpoint_psi = c10.psi( 8 );
		Real const right_endpoint_omega = c10.omega( 8 );

		// attempt a SegmentInsert
		SegmentInsert si( Interval( 4, 8 ), "EEEE^E", "GGGG^G", h10, true, N );
		si.modify( c10 );

		// force coordinate update
		c10.residue( 7 );

		// check post-conditions
		TS_ASSERT_EQUALS( si.flanking_left_nres(), 4 );
		TS_ASSERT_EQUALS( si.flanking_right_nres(), 1 );
		TS_ASSERT_EQUALS( si.flanking_left_ss(), "EEEE" );
		TS_ASSERT_EQUALS( si.flanking_right_ss(), "E" );
		TS_ASSERT_EQUALS( si.flanking_left_aa(), "GGGG" );
		TS_ASSERT_EQUALS( si.flanking_right_aa(), "G" );

		TS_ASSERT_EQUALS( si.interval().left, 4 );
		TS_ASSERT_EQUALS( si.interval().right, 18 );
		TS_ASSERT( !si.original_interval_valid() );

		TS_ASSERT_EQUALS( c10.size(), 20 );
		TS_ASSERT_EQUALS( c10.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( c10.annotated_sequence(), "R[ARG:NtermProteinFull]DDGGGGKEEEEEEEEKGDR[ARG:CtermProteinFull]" );
		TS_ASSERT_EQUALS( c10.secstruct(), "LLLEEEEHHHHHHHHHHELL" );

		TS_ASSERT_EQUALS( c10.residue( 3 ).xyz( "CA" ), left_minus_one_CA );
		TS_ASSERT_EQUALS( c10.residue( 19 ).xyz( "CA" ), right_plus_one_CA );

		for ( Size i = 8; i <= 17; ++i ) { // check insert phi
			// 8 is the first residue of the insert so the phi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 8 ) {
				TS_ASSERT_DELTA( c10.phi( i ), -60.0, 0.001 );
			}

			// 17 is the last residue of the insert so the psi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 17 ) {
				TS_ASSERT_DELTA( c10.psi( i ), -45.0, 0.001 );
			}
		}

		for ( Size i = 4; i <= 18; ++i ) { // check flanking + insert omega
			if ( !c10.fold_tree().is_cutpoint( i ) ) {
				TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( c10.omega( i ) ), 180.0, 0.001 );
			}
		}

		// check bb torsions at endpoints, we requested that SegmentInsert retain
		// them
		TS_ASSERT_DELTA( c10.phi( si.interval().left ), left_endpoint_phi, 0.001 );
		TS_ASSERT_DELTA( c10.psi( si.interval().right ), right_endpoint_psi, 0.001 );
		TS_ASSERT_DELTA( c10.omega( si.interval().right ), right_endpoint_omega, 0.001 );

		// N-side connection, so make sure a cutpoint exists somewhere on the right
		TS_ASSERT( c10.fold_tree().is_cutpoint( 17 ) ); // 17 should be a cutpoint

		TS_ASSERT_DELTA( c10.residue( 8 ).xyz( "CA" ).distance( c10.residue( 12 ).xyz( "CA" ) ), f_m_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 12 ).xyz( "CA" ).distance( c10.residue( 17 ).xyz( "CA" ) ), m_l_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 8 ).xyz( "CA" ).distance( c10.residue( 17 ).xyz( "CA" ) ), f_l_dist, 0.001 );
	}


	/// @brief test insertion with C-side connection, simple replacement, completely
	///  internal insertion
	void test_C_side_connection_simple_internal() {
		using namespace protocols::forge::build::SegmentInsertConnectionScheme;

		using numeric::principal_angle_degrees;
		using numeric::nonnegative_principal_angle_degrees;

		// create dummy poses
		Pose h10 = helix10_distorted_pose();
		Pose c10 = continuous10_pose();

		// Store CA <-> CA distances between first, middle, and last residues
		// of h10 pairwise.  These shouldn't change after the operation.
		Real const f_m_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 5 ).xyz( "CA" ) );
		Real const m_l_dist = h10.residue( 5 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );
		Real const f_l_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );

		// Store CA positions of flanking left-1 and right+1 residues.
		// These shouldn't change after the operation.
		Vector const left_minus_one_CA = c10.residue( 3 ).xyz( "CA" );
		Vector const right_plus_one_CA = c10.residue( 9 ).xyz( "CA" );

		// Store bb torsions at endpoints; we are requesting that SegmentInsert
		// retain them.
		Real const left_endpoint_phi = c10.phi( 4 );
		Real const right_endpoint_psi = c10.psi( 8 );
		Real const right_endpoint_omega = c10.omega( 8 );

		// attempt a SegmentInsert
		SegmentInsert si( Interval( 4, 8 ), "EEEE^E", "GGGG^G", h10, true, C );
		si.modify( c10 );

		// force coordinate update
		c10.residue( 7 );

		// check post-conditions
		TS_ASSERT_EQUALS( si.flanking_left_nres(), 4 );
		TS_ASSERT_EQUALS( si.flanking_right_nres(), 1 );
		TS_ASSERT_EQUALS( si.flanking_left_ss(), "EEEE" );
		TS_ASSERT_EQUALS( si.flanking_right_ss(), "E" );
		TS_ASSERT_EQUALS( si.flanking_left_aa(), "GGGG" );
		TS_ASSERT_EQUALS( si.flanking_right_aa(), "G" );

		TS_ASSERT_EQUALS( si.interval().left, 4 );
		TS_ASSERT_EQUALS( si.interval().right, 18 );
		TS_ASSERT( !si.original_interval_valid() );

		TS_ASSERT_EQUALS( c10.size(), 20 );
		TS_ASSERT_EQUALS( c10.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( c10.annotated_sequence(), "R[ARG:NtermProteinFull]DDGGGGKEEEEEEEEKGDR[ARG:CtermProteinFull]" );
		TS_ASSERT_EQUALS( c10.secstruct(), "LLLEEEEHHHHHHHHHHELL" );

		TS_ASSERT_EQUALS( c10.residue( 3 ).xyz( "CA" ), left_minus_one_CA );
		TS_ASSERT_EQUALS( c10.residue( 19 ).xyz( "CA" ), right_plus_one_CA );

		for ( Size i = 8; i <= 17; ++i ) { // check insert phi
			// 8 is the first residue of the insert so the phi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 8 ) {
				TS_ASSERT_DELTA( c10.phi( i ), -60.0, 0.001 );
			}

			// 17 is the last residue of the insert so the psi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 17 ) {
				TS_ASSERT_DELTA( c10.psi( i ), -45.0, 0.001 );
			}
		}

		for ( Size i = 4; i <= 18; ++i ) { // check flanking + insert omega
			if ( !c10.fold_tree().is_cutpoint( i ) ) {
				TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( c10.omega( i ) ), 180.0, 0.001 );
			}
		}

		// check bb torsions at endpoints, we requested that SegmentInsert retain
		// them
		TS_ASSERT_DELTA( c10.phi( si.interval().left ), left_endpoint_phi, 0.001 );
		TS_ASSERT_DELTA( c10.psi( si.interval().right ), right_endpoint_psi, 0.001 );
		TS_ASSERT_DELTA( c10.omega( si.interval().right ), right_endpoint_omega, 0.001 );

		// make sure a single cutpoint exists somewhere on the left
		Size n_left_cut = 0;
		for ( Size i = 3; i <= 7; ++i ) {
			if ( c10.fold_tree().is_cutpoint( i ) ) {
				++n_left_cut;
			}
		}
		TS_ASSERT_EQUALS( n_left_cut, 1 );

		TS_ASSERT_DELTA( c10.residue( 8 ).xyz( "CA" ).distance( c10.residue( 12 ).xyz( "CA" ) ), f_m_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 12 ).xyz( "CA" ).distance( c10.residue( 17 ).xyz( "CA" ) ), m_l_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 8 ).xyz( "CA" ).distance( c10.residue( 17 ).xyz( "CA" ) ), f_l_dist, 0.001 );
	}


	/// @brief test insertion with N-side connection, simple replacement, completely
	///  internal insertion, pure insertion
	void test_N_side_connection_simple_internal_pure() {
		using namespace protocols::forge::build::SegmentInsertConnectionScheme;

		using numeric::principal_angle_degrees;
		using numeric::nonnegative_principal_angle_degrees;

		// create dummy poses
		Pose h10 = helix10_distorted_pose();
		Pose c10 = continuous10_pose();

		// Store CA <-> CA distances between first, middle, and last residues
		// of h10 pairwise.  These shouldn't change after the operation.
		Real const f_m_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 5 ).xyz( "CA" ) );
		Real const m_l_dist = h10.residue( 5 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );
		Real const f_l_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );

		// Store CA positions of flanking left-1 and right+1 residues.
		// These shouldn't change after the operation.
		Vector const left_minus_one_CA = c10.residue( 7 ).xyz( "CA" );
		Vector const right_plus_one_CA = c10.residue( 8 ).xyz( "CA" );

		// attempt a SegmentInsert
		SegmentInsert si( Interval( 0, 7 ), "EEEE^E", "GGGG^G", h10, false, N );
		si.modify( c10 );

		// force coordinate update
		c10.residue( 7 );

		// check post-conditions
		TS_ASSERT_EQUALS( si.flanking_left_nres(), 4 );
		TS_ASSERT_EQUALS( si.flanking_right_nres(), 1 );
		TS_ASSERT_EQUALS( si.flanking_left_ss(), "EEEE" );
		TS_ASSERT_EQUALS( si.flanking_right_ss(), "E" );
		TS_ASSERT_EQUALS( si.flanking_left_aa(), "GGGG" );
		TS_ASSERT_EQUALS( si.flanking_right_aa(), "G" );

		TS_ASSERT_EQUALS( si.interval().left, 8 );
		TS_ASSERT_EQUALS( si.interval().right, 22 );
		TS_ASSERT( !si.original_interval_valid() );

		TS_ASSERT_EQUALS( c10.size(), 25 );
		TS_ASSERT_EQUALS( c10.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( c10.annotated_sequence(), "R[ARG:NtermProteinFull]DDDDDDGGGGKEEEEEEEEKGDDR[ARG:CtermProteinFull]" );
		TS_ASSERT_EQUALS( c10.secstruct(), "LLLLLLLEEEEHHHHHHHHHHELLL" );

		TS_ASSERT_EQUALS( c10.residue( 7 ).xyz( "CA" ), left_minus_one_CA );
		TS_ASSERT_EQUALS( c10.residue( 23 ).xyz( "CA" ), right_plus_one_CA );

		for ( Size i = 12; i <= 21; ++i ) { // check insert phi
			// 12 is the first residue of the insert so the phi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 12 ) {
				TS_ASSERT_DELTA( c10.phi( i ), -60.0, 0.001 );
			}

			// 21 is the last residue of the insert so the psi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 21 ) {
				TS_ASSERT_DELTA( c10.psi( i ), -45.0, 0.001 );
			}
		}

		for ( Size i = 8; i <= 22; ++i ) { // check flanking + insert omega
			if ( !c10.fold_tree().is_cutpoint( i ) ) {
				TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( c10.omega( i ) ), 180.0, 0.001 );
			}
		}

		// N-side connection, so make sure a cutpoint exists somewhere on the right
		Size n_right_cut = 0;
		for ( Size i = 21; i <= 22; ++i ) {
			if ( c10.fold_tree().is_cutpoint( i ) ) {
				++n_right_cut;
			}
		}
		TS_ASSERT_EQUALS( n_right_cut, 1 );

		TS_ASSERT_DELTA( c10.residue( 12 ).xyz( "CA" ).distance( c10.residue( 16 ).xyz( "CA" ) ), f_m_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 16 ).xyz( "CA" ).distance( c10.residue( 21 ).xyz( "CA" ) ), m_l_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 12 ).xyz( "CA" ).distance( c10.residue( 21 ).xyz( "CA" ) ), f_l_dist, 0.001 );
	}


	/// @brief test insertion with C-side connection, simple replacement, completely
	///  internal insertion, pure insertion
	void test_C_side_connection_simple_internal_pure() {
		using namespace protocols::forge::build::SegmentInsertConnectionScheme;

		using numeric::principal_angle_degrees;
		using numeric::nonnegative_principal_angle_degrees;

		// create dummy poses
		Pose h10 = helix10_distorted_pose();
		Pose c10 = continuous10_pose();

		// Store CA <-> CA distances between first, middle, and last residues
		// of h10 pairwise.  These shouldn't change after the operation.
		Real const f_m_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 5 ).xyz( "CA" ) );
		Real const m_l_dist = h10.residue( 5 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );
		Real const f_l_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );

		// Store CA positions of flanking left-1 and right+1 residues.
		// These shouldn't change after the operation.
		Vector const left_minus_one_CA = c10.residue( 7 ).xyz( "CA" );
		Vector const right_plus_one_CA = c10.residue( 8 ).xyz( "CA" );

		// attempt a SegmentInsert
		SegmentInsert si( Interval( 0, 7 ), "EEEE^E", "GGGG^G", h10, false, C );
		si.modify( c10 );

		// force coordinate update
		c10.residue( 7 );

		// check post-conditions
		TS_ASSERT_EQUALS( si.flanking_left_nres(), 4 );
		TS_ASSERT_EQUALS( si.flanking_right_nres(), 1 );
		TS_ASSERT_EQUALS( si.flanking_left_ss(), "EEEE" );
		TS_ASSERT_EQUALS( si.flanking_right_ss(), "E" );
		TS_ASSERT_EQUALS( si.flanking_left_aa(), "GGGG" );
		TS_ASSERT_EQUALS( si.flanking_right_aa(), "G" );

		TS_ASSERT_EQUALS( si.interval().left, 8 );
		TS_ASSERT_EQUALS( si.interval().right, 22 );
		TS_ASSERT( !si.original_interval_valid() );

		TS_ASSERT_EQUALS( c10.size(), 25 );
		TS_ASSERT_EQUALS( c10.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( c10.annotated_sequence(), "R[ARG:NtermProteinFull]DDDDDDGGGGKEEEEEEEEKGDDR[ARG:CtermProteinFull]" );
		TS_ASSERT_EQUALS( c10.secstruct(), "LLLLLLLEEEEHHHHHHHHHHELLL" );

		TS_ASSERT_EQUALS( c10.residue( 7 ).xyz( "CA" ), left_minus_one_CA );
		TS_ASSERT_EQUALS( c10.residue( 23 ).xyz( "CA" ), right_plus_one_CA );

		for ( Size i = 12; i <= 21; ++i ) { // check insert phi
			// 12 is the first residue of the insert so the phi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 12 ) {
				TS_ASSERT_DELTA( c10.phi( i ), -60.0, 0.001 );
			}

			// 21 is the last residue of the insert so the psi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 21 ) {
				TS_ASSERT_DELTA( c10.psi( i ), -45.0, 0.001 );
			}
		}

		for ( Size i = 8; i <= 22; ++i ) { // check flanking + insert omega
			if ( !c10.fold_tree().is_cutpoint( i ) ) {
				TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( c10.omega( i ) ), 180.0, 0.001 );
			}
		}

		// N-side connection, so make sure a cutpoint exists somewhere on the left
		Size n_left_cut = 0;
		for ( Size i = 7; i <= 11; ++i ) {
			if ( c10.fold_tree().is_cutpoint( i ) ) {
				++n_left_cut;
			}
		}
		TS_ASSERT_EQUALS( n_left_cut, 1 );

		TS_ASSERT_DELTA( c10.residue( 12 ).xyz( "CA" ).distance( c10.residue( 16 ).xyz( "CA" ) ), f_m_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 16 ).xyz( "CA" ).distance( c10.residue( 21 ).xyz( "CA" ) ), m_l_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 12 ).xyz( "CA" ).distance( c10.residue( 21 ).xyz( "CA" ) ), f_l_dist, 0.001 );
	}


	/// @brief test insertion with N-side connection, simple replacement, n-terminal insertion
	void test_C_side_connection_simple_n_term() {
		using namespace protocols::forge::build::SegmentInsertConnectionScheme;

		using numeric::principal_angle_degrees;
		using numeric::nonnegative_principal_angle_degrees;

		// create dummy poses
		Pose h10 = helix10_distorted_pose();
		Pose c10 = continuous10_pose();

		// Store CA <-> CA distances between first, middle, and last residues
		// of h10 pairwise.  These shouldn't change after the operation.
		Real const f_m_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 5 ).xyz( "CA" ) );
		Real const m_l_dist = h10.residue( 5 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );
		Real const f_l_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );

		// Store CA positions of flanking right+1 residue.
		// This shouldn't change after the operation.
		Vector const right_plus_one_CA = c10.residue( 1 ).xyz( "CA" );

		// attempt a SegmentInsert
		SegmentInsert si( Interval( 0, 0 ), "EEEE^E", "GGGG^G", h10, false, N );
		si.modify( c10 );

		// force coordinate update
		c10.residue( 7 );

		// check post-conditions
		TS_ASSERT_EQUALS( si.flanking_left_nres(), 4 );
		TS_ASSERT_EQUALS( si.flanking_right_nres(), 1 );
		TS_ASSERT_EQUALS( si.flanking_left_ss(), "EEEE" );
		TS_ASSERT_EQUALS( si.flanking_right_ss(), "E" );
		TS_ASSERT_EQUALS( si.flanking_left_aa(), "GGGG" );
		TS_ASSERT_EQUALS( si.flanking_right_aa(), "G" );

		TS_ASSERT_EQUALS( si.interval().left, 1 );
		TS_ASSERT_EQUALS( si.interval().right, 15 );
		TS_ASSERT( !si.original_interval_valid() );

		TS_ASSERT_EQUALS( c10.size(), 25 );
		TS_ASSERT_EQUALS( c10.fold_tree().num_cutpoint(), 0 );
		TS_ASSERT_EQUALS( c10.annotated_sequence(), "G[GLY:NtermProteinFull]GGGKEEEEEEEEKGRDDDDDDDDR[ARG:CtermProteinFull]" );
		TS_ASSERT_EQUALS( c10.secstruct(), "EEEEHHHHHHHHHHELLLLLLLLLL" );

		TS_ASSERT_EQUALS( c10.residue( 16 ).xyz( "CA" ), right_plus_one_CA );

		for ( Size i = 5; i <= 14; ++i ) { // check insert phi
			// 5 is the first residue of the insert so the phi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 5 ) {
				TS_ASSERT_DELTA( c10.phi( i ), -60.0, 0.001 );
			}

			// 14 is the last residue of the insert so the psi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 14 ) {
				TS_ASSERT_DELTA( c10.psi( i ), -45.0, 0.001 );
			}
		}

		for ( Size i = 1; i <= 15; ++i ) { // check flanking + insert omega
			if ( !c10.fold_tree().is_cutpoint( i ) ) {
				TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( c10.omega( i ) ), 180.0, 0.001 );
			}
		}

		TS_ASSERT_DELTA( c10.residue( 5 ).xyz( "CA" ).distance( c10.residue( 9 ).xyz( "CA" ) ), f_m_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 9 ).xyz( "CA" ).distance( c10.residue( 14 ).xyz( "CA" ) ), m_l_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 5 ).xyz( "CA" ).distance( c10.residue( 14 ).xyz( "CA" ) ), f_l_dist, 0.001 );
	}


	/// @brief test insertion with N-side connection, simple replacement, c-terminal insertion
	void test_N_side_connection_simple_c_term() {
		using namespace protocols::forge::build::SegmentInsertConnectionScheme;

		using numeric::principal_angle_degrees;
		using numeric::nonnegative_principal_angle_degrees;

		// create dummy poses
		Pose h10 = helix10_distorted_pose();
		Pose c10 = continuous10_pose();

		// Store CA <-> CA distances between first, middle, and last residues
		// of h10 pairwise.  These shouldn't change after the operation.
		Real const f_m_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 5 ).xyz( "CA" ) );
		Real const m_l_dist = h10.residue( 5 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );
		Real const f_l_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );

		// Store CA positions of flanking left-1 residue.
		// This shouldn't change after the operation.
		Vector const left_minus_one_CA = c10.residue( 10 ).xyz( "CA" );

		// attempt a SegmentInsert
		SegmentInsert si( Interval( 0, c10.size() ), "E^EEEE", "G^GGGG", h10, false, C );
		si.modify( c10 );

		// force coordinate update
		c10.residue( 7 );

		// check post-conditions
		TS_ASSERT_EQUALS( si.flanking_left_nres(), 1 );
		TS_ASSERT_EQUALS( si.flanking_right_nres(), 4 );
		TS_ASSERT_EQUALS( si.flanking_left_ss(), "E" );
		TS_ASSERT_EQUALS( si.flanking_right_ss(), "EEEE" );
		TS_ASSERT_EQUALS( si.flanking_left_aa(), "G" );
		TS_ASSERT_EQUALS( si.flanking_right_aa(), "GGGG" );

		TS_ASSERT_EQUALS( si.interval().left, 11 );
		TS_ASSERT_EQUALS( si.interval().right, 25 );
		TS_ASSERT( !si.original_interval_valid() );

		TS_ASSERT_EQUALS( c10.size(), 25 );
		TS_ASSERT_EQUALS( c10.fold_tree().num_cutpoint(), 0 );
		TS_ASSERT_EQUALS( c10.annotated_sequence(), "R[ARG:NtermProteinFull]DDDDDDDDRGKEEEEEEEEKGGGG[GLY:CtermProteinFull]" );
		TS_ASSERT_EQUALS( c10.secstruct(), "LLLLLLLLLLEHHHHHHHHHHEEEE" );

		TS_ASSERT_EQUALS( c10.residue( 10 ).xyz( "CA" ), left_minus_one_CA );

		for ( Size i = 12; i <= 21; ++i ) { // check insert phi
			// 12 is the first residue of the insert so the phi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 12 ) {
				TS_ASSERT_DELTA( c10.phi( i ), -60.0, 0.001 );
			}

			// 21 is the last residue of the insert so the psi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 21 ) {
				TS_ASSERT_DELTA( c10.psi( i ), -45.0, 0.001 );
			}
		}

		for ( Size i = 11; i <= 21; ++i ) { // check flanking + insert omega
			if ( !c10.fold_tree().is_cutpoint( i ) ) {
				TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( c10.omega( i ) ), 180.0, 0.001 );
			}
		}

		TS_ASSERT_DELTA( c10.residue( 12 ).xyz( "CA" ).distance( c10.residue( 16 ).xyz( "CA" ) ), f_m_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 16 ).xyz( "CA" ).distance( c10.residue( 21 ).xyz( "CA" ) ), m_l_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 12 ).xyz( "CA" ).distance( c10.residue( 21 ).xyz( "CA" ) ), f_l_dist, 0.001 );
	}


	/// @brief test insertion with N-side connection, simple replacement,
	///  completely internal, but no flanking region on the left
	void test_N_side_connection_no_left_flanking_region() {
		using namespace protocols::forge::build::SegmentInsertConnectionScheme;

		using numeric::principal_angle_degrees;
		using numeric::nonnegative_principal_angle_degrees;

		// create dummy poses
		Pose h10 = helix10_distorted_pose();
		Pose c10 = continuous10_pose();

		// Store CA <-> CA distances between first, middle, and last residues
		// of h10 pairwise.  These shouldn't change after the operation.
		Real const f_m_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 5 ).xyz( "CA" ) );
		Real const m_l_dist = h10.residue( 5 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );
		Real const f_l_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );

		// Store CA positions of flanking left-1 and right+1 residues.
		// These shouldn't change after the operation.
		Vector const left_minus_one_CA = c10.residue( 3 ).xyz( "CA" );
		Vector const right_plus_one_CA = c10.residue( 9 ).xyz( "CA" );

		// Store bb torsions at endpoints; we are requesting that SegmentInsert
		// retain them.
		Real const left_endpoint_phi = c10.phi( 4 );
		Real const right_endpoint_psi = c10.psi( 8 );
		Real const right_endpoint_omega = c10.omega( 8 );

		// attempt a SegmentInsert
		SegmentInsert si( Interval( 4, 8 ), "^EEEEE", "^GGGGG", h10, true, N );
		si.modify( c10 );

		// force coordinate update
		c10.residue( 7 );

		// check post-conditions
		TS_ASSERT_EQUALS( si.flanking_left_nres(), 0 );
		TS_ASSERT_EQUALS( si.flanking_right_nres(), 5 );
		TS_ASSERT_EQUALS( si.flanking_left_ss(), "" );
		TS_ASSERT_EQUALS( si.flanking_right_ss(), "EEEEE" );
		TS_ASSERT_EQUALS( si.flanking_left_aa(), "" );
		TS_ASSERT_EQUALS( si.flanking_right_aa(), "GGGGG" );

		TS_ASSERT_EQUALS( si.interval().left, 4 );
		TS_ASSERT_EQUALS( si.interval().right, 18 );
		TS_ASSERT( !si.original_interval_valid() );

		TS_ASSERT_EQUALS( c10.size(), 20 );
		TS_ASSERT_EQUALS( c10.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( c10.annotated_sequence(), "R[ARG:NtermProteinFull]DDKEEEEEEEEKGGGGGDR[ARG:CtermProteinFull]" );
		TS_ASSERT_EQUALS( c10.secstruct(), "LLLHHHHHHHHHHEEEEELL" );

		TS_ASSERT_EQUALS( c10.residue( 3 ).xyz( "CA" ), left_minus_one_CA );
		TS_ASSERT_EQUALS( c10.residue( 19 ).xyz( "CA" ), right_plus_one_CA );

		for ( Size i = 4; i <= 13; ++i ) { // check insert phi
			// 4 is the first residue of the insert so the psi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 4 ) {
				TS_ASSERT_DELTA( c10.phi( i ), -60.0, 0.001 );
			}

			// 17 is the last residue of the insert so the psi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 13 ) {
				TS_ASSERT_DELTA( c10.psi( i ), -45.0, 0.001 );
			}
		}

		for ( Size i = 4; i <= 18; ++i ) { // check flanking + insert omega
			if ( !c10.fold_tree().is_cutpoint( i ) ) {
				TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( c10.omega( i ) ), 180.0, 0.001 );
			}
		}

		// check bb torsions at endpoints, we requested that SegmentInsert retain
		// them
		TS_ASSERT_DELTA( c10.phi( si.interval().left ), left_endpoint_phi, 0.001 );
		TS_ASSERT_DELTA( c10.psi( si.interval().right ), right_endpoint_psi, 0.001 );
		TS_ASSERT_DELTA( c10.omega( si.interval().right ), right_endpoint_omega, 0.001 );

		// N-side connection, so make sure a cutpoint exists somewhere on the right
		bool found_cutpoint = false;
		for ( Size i = 13; i <= 18 && !found_cutpoint; ++i ) {
			if ( c10.fold_tree().is_cutpoint( i ) ) {
				found_cutpoint = true;
			}
		}
		TS_ASSERT( found_cutpoint );

		TS_ASSERT_DELTA( c10.residue( 4 ).xyz( "CA" ).distance( c10.residue( 8 ).xyz( "CA" ) ), f_m_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 8 ).xyz( "CA" ).distance( c10.residue( 13 ).xyz( "CA" ) ), m_l_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 4 ).xyz( "CA" ).distance( c10.residue( 13 ).xyz( "CA" ) ), f_l_dist, 0.001 );
	}


	/// @brief test insertion with C-side connection, simple replacement,
	///  completely internal, but no flanking region on the right
	void test_N_side_connection_no_right_flanking_region() {
		using namespace protocols::forge::build::SegmentInsertConnectionScheme;

		using numeric::principal_angle_degrees;
		using numeric::nonnegative_principal_angle_degrees;

		// create dummy poses
		Pose h10 = helix10_distorted_pose();
		Pose c10 = continuous10_pose();

		// Store CA <-> CA distances between first, middle, and last residues
		// of h10 pairwise.  These shouldn't change after the operation.
		Real const f_m_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 5 ).xyz( "CA" ) );
		Real const m_l_dist = h10.residue( 5 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );
		Real const f_l_dist = h10.residue( 1 ).xyz( "CA" ).distance( h10.residue( 10 ).xyz( "CA" ) );

		// Store CA positions of flanking left-1 and right+1 residues.
		// These shouldn't change after the operation.
		Vector const left_minus_one_CA = c10.residue( 3 ).xyz( "CA" );
		Vector const right_plus_one_CA = c10.residue( 9 ).xyz( "CA" );

		// Store bb torsions at endpoints; we are requesting that SegmentInsert
		// retain them.
		Real const left_endpoint_phi = c10.phi( 4 );
		Real const right_endpoint_psi = c10.psi( 8 );
		Real const right_endpoint_omega = c10.omega( 8 );

		// attempt a SegmentInsert
		SegmentInsert si( Interval( 4, 8 ), "EEEEE^", "GGGGG^", h10, true, C );
		si.modify( c10 );

		// force coordinate update
		c10.residue( 7 );

		// check post-conditions
		TS_ASSERT_EQUALS( si.flanking_left_nres(), 5 );
		TS_ASSERT_EQUALS( si.flanking_right_nres(), 0 );
		TS_ASSERT_EQUALS( si.flanking_left_ss(), "EEEEE" );
		TS_ASSERT_EQUALS( si.flanking_right_ss(), "" );
		TS_ASSERT_EQUALS( si.flanking_left_aa(), "GGGGG" );
		TS_ASSERT_EQUALS( si.flanking_right_aa(), "" );

		TS_ASSERT_EQUALS( si.interval().left, 4 );
		TS_ASSERT_EQUALS( si.interval().right, 18 );
		TS_ASSERT( !si.original_interval_valid() );

		TS_ASSERT_EQUALS( c10.size(), 20 );
		TS_ASSERT_EQUALS( c10.fold_tree().num_cutpoint(), 1 );
		TS_ASSERT_EQUALS( c10.annotated_sequence(), "R[ARG:NtermProteinFull]DDGGGGGKEEEEEEEEKDR[ARG:CtermProteinFull]" );
		TS_ASSERT_EQUALS( c10.secstruct(), "LLLEEEEEHHHHHHHHHHLL" );

		TS_ASSERT_EQUALS( c10.residue( 3 ).xyz( "CA" ), left_minus_one_CA );
		TS_ASSERT_EQUALS( c10.residue( 19 ).xyz( "CA" ), right_plus_one_CA );

		for ( Size i = 9; i <= 18; ++i ) { // check insert phi
			// 9 is the first residue of the insert so the phi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 9 ) {
				TS_ASSERT_DELTA( c10.phi( i ), -60.0, 0.001 );
			}

			// 18 is the last residue of the insert so the psi for this doesn't
			// really exist in the original h10 pose; skip it
			if ( i != 18 ) {
				TS_ASSERT_DELTA( c10.psi( i ), -45.0, 0.001 );
			}
		}

		for ( Size i = 4; i <= 18; ++i ) { // check flanking + insert omega
			if ( !c10.fold_tree().is_cutpoint( i ) ) {
				TS_ASSERT_DELTA( nonnegative_principal_angle_degrees( c10.omega( i ) ), 180.0, 0.001 );
			}
		}

		// check bb torsions at endpoints, we requested that SegmentInsert retain
		// them
		TS_ASSERT_DELTA( c10.phi( si.interval().left ), left_endpoint_phi, 0.001 );
		TS_ASSERT_DELTA( c10.psi( si.interval().right ), right_endpoint_psi, 0.001 );
		TS_ASSERT_DELTA( c10.omega( si.interval().right ), right_endpoint_omega, 0.001 );

		// make sure a single cutpoint exists somewhere on the left
		Size n_left_cut = 0;
		for ( Size i = 3; i <= 8; ++i ) {
			if ( c10.fold_tree().is_cutpoint( i ) ) {
				++n_left_cut;
			}
		}
		TS_ASSERT_EQUALS( n_left_cut, 1 );

		TS_ASSERT_DELTA( c10.residue( 9 ).xyz( "CA" ).distance( c10.residue( 13 ).xyz( "CA" ) ), f_m_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 13 ).xyz( "CA" ).distance( c10.residue( 18 ).xyz( "CA" ) ), m_l_dist, 0.001 );
		TS_ASSERT_DELTA( c10.residue( 9 ).xyz( "CA" ).distance( c10.residue( 18 ).xyz( "CA" ) ), f_l_dist, 0.001 );
	}


};
