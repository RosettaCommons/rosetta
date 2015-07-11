// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/FoldTree.cxxtest.hh
/// @brief  test suite for core::kinematics::FoldTree.cc
/// @author Christopher Miles (cmiles@uw.edu)
/// @author Steven Lewis smlewi@gmail.com test_utility_function_remodel_fold_tree_to_account_for_insertion

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>

// Package Header
#include <core/kinematics/util.hh>  //for test_utility_function_remodel_fold_tree_to_account_for_insertion

// Project headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

//Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <sstream>

using core::kinematics::FoldTree;
using core::pose::Pose;

class FoldTreeTest : public CxxTest::TestSuite {
public:
	Pose pose_;

	void setUp() {
		core_init();
		core::import_pose::pose_from_pdb( pose_, "core/kinematics/test.pdb" );
	}

	// Ensure that repeated calls to hash value yield identical results when applied to the same FoldTree.
	void test_hash_value_unmodified() {
		TS_ASSERT_EQUALS( pose_.fold_tree().hash_value(), pose_.fold_tree().hash_value() );
	}

	// Ensure that operations that modify the FoldTree have an effect on the calculation of the hash value.
	void test_hash_value_modified() {
		FoldTree mod_tree( pose_.fold_tree() );
		mod_tree.new_jump( 1, 3, 2 );

		size_t hash_orig = pose_.fold_tree().hash_value();
		size_t hash_mod = mod_tree.hash_value();
		TS_ASSERT_DIFFERS( hash_orig, hash_mod );
	}

	// A quick check to see whether the FoldTree considers the final residues in a simple fold tree to be a cutpoint.
	void test_last_residue_is_cutpoint() {
		FoldTree const & tree( pose_.fold_tree() );
		TS_ASSERT( tree.is_cutpoint( pose_.total_residue() ) );
	}

	void test_boundary_right() {
		FoldTree const & tree = pose_.fold_tree();
		unsigned const n = pose_.total_residue();

		for ( unsigned i = 1; i <= n; ++i ) {
			if ( tree.is_root( i ) ) { continue; }

			// Latter is a shortcut to the former
			TS_ASSERT_EQUALS( n, tree.get_residue_edge( i ).stop() );
			TS_ASSERT_EQUALS( n, tree.boundary_right( i ) );
		}
	}

	void test_boundary_left() {
		FoldTree const & tree = pose_.fold_tree();

		for ( unsigned i = 1; i < pose_.total_residue(); ++i ) {
			if ( tree.is_root( i ) || i == tree.nres() ) { continue; }

			// Latter is a shortcut to the former
			TS_ASSERT_EQUALS( 1, tree.get_residue_edge( i ).start() );
			TS_ASSERT_EQUALS( 1, tree.boundary_left( i ) );
		}
	}

	void test_get_chemical_edges() {
		FoldTree ft;
		ft.add_edge( 1, 5, -1 );
		ft.add_edge( 3, 6, "ND2", "C1" );
		ft.add_edge( 6, 10, -1 );
		ft.add_edge( 1, 11, 1 );
		ft.add_edge( 11, 20, -1 );
		ft.add_edge( 15, 21, "O2", "C1" );
		ft.add_edge( 21, 23, -1 );
		ft.add_edge( 20, 24, 2 );
		ft.add_edge( 24, 30, -1 );
		ft.add_edge( 22, 31, "O6", "C1" );
		ft.add_edge( 31, 33, -1 );

		TS_ASSERT( ft.check_fold_tree() );

		utility::vector1< core::kinematics::Edge > const chemical_edges( ft.get_chemical_edges() );
		TS_ASSERT_EQUALS( chemical_edges.size(), 3 );

		core::kinematics::Edge const edge1( chemical_edges[ 1 ] );
		core::kinematics::Edge const edge2( chemical_edges[ 2 ] );
		core::kinematics::Edge const edge3( chemical_edges[ 3 ] );

		TS_ASSERT( edge1.valid() );
		TS_ASSERT( edge1.is_chemical_bond() );
		TS_ASSERT( ! edge1.is_jump() );  //assert that it is not a jump
		TS_ASSERT( ! edge1.is_peptide() );  //assert that it is not a peptide

		TS_ASSERT_EQUALS( edge2.start_atom(), "O2" );
		TS_ASSERT_EQUALS( edge2.stop_atom(), "C1" );
		TS_ASSERT_EQUALS( edge2.label(), -2 );

		TS_ASSERT_EQUALS( edge3.start(), 22 );
		TS_ASSERT_EQUALS( edge3.stop(), 31 );
	}


	// SML 10-19-12 this is a utility function in kinematics/util.*, but it seems appropriate to test with FoldTree
	void test_utility_function_remodel_fold_tree_to_account_for_insertion() {
		std::istringstream reusable_FT_istream;

		// Test a simple case.
		FoldTree const fifty_foldtree( 50 ), sixty_foldtree( 60 );
		// should have one edge, 1-50
		FoldTree fifty_foldtree_after(
				core::kinematics::remodel_fold_tree_to_account_for_insertion( fifty_foldtree, 25, 10 ) );
		// should have one edge, 1-60
		TS_ASSERT_EQUALS( fifty_foldtree_after, sixty_foldtree );

		// fifty_foldtree FOLD_TREE  EDGE 1 50 -1
		// sixty_foldtree FOLD_TREE  EDGE 1 60 -1
		// fifty_foldtree_after FOLD_TREE  EDGE 1 60 -1
		// *************************************************C
		// ***********************************************************C
		// ***********************************************************C

		//Test a more complex case:
		using core::kinematics::Edge;
		core::Size const cbreak( 25 );  //arbitrary chainbreak position
		FoldTree twochain_tree( 50 );
		twochain_tree.clear();
		twochain_tree.add_edge( Edge( 1, cbreak, Edge::PEPTIDE ) );
		twochain_tree.add_edge( Edge( cbreak + 1, 50, Edge::PEPTIDE ) );
		twochain_tree.add_edge( Edge( cbreak, cbreak + 1, 1 ) );
		twochain_tree.reorder( 1 );

		// twochain_tree FOLD_TREE  EDGE 1 25 -1  EDGE 25 26 1  EDGE 26 50 -1
		// ************************1/C1***********************C

		FoldTree add_before_jump;
		reusable_FT_istream.clear();
		reusable_FT_istream.str( "FOLD_TREE  EDGE 1 35 -1  EDGE 35 36 1  EDGE 36 60 -1" );
		reusable_FT_istream >> add_before_jump;

		TS_ASSERT_EQUALS(add_before_jump,
				core::kinematics::remodel_fold_tree_to_account_for_insertion( twochain_tree, 20, 10 ) );

		core::kinematics::FoldTree add_after_jump;
		reusable_FT_istream.clear();
		reusable_FT_istream.str("FOLD_TREE  EDGE 1 25 -1  EDGE 25 26 1  EDGE 26 60 -1");
		reusable_FT_istream >> add_after_jump;

		TS_ASSERT_EQUALS(add_after_jump,
			core::kinematics::remodel_fold_tree_to_account_for_insertion(twochain_tree, 30, 10));


		// This next line tests that it crashes (as it should) if you try to insert on a jumping point using this
		// simple algorithm
		TS_ASSERT_THROWS_ANYTHING(
				core::kinematics::remodel_fold_tree_to_account_for_insertion( twochain_tree, 25, 10 ) );


		//Test complex case #1 (stolen from AnchoredDesign integration test)
		core::kinematics::FoldTree anchored_design_fold_tree;
		reusable_FT_istream.clear();
		reusable_FT_istream.str(
				"FOLD_TREE EDGE 1 59 -1  JEDGE 59 85 1 C N  INTRA_RES_STUB  EDGE 85 78 -1  EDGE 85 90 -1  "
				"EDGE 90 132 -1  EDGE 132 133 -1  EDGE 90 76 2  EDGE 132 143 3  EDGE 76 60 -1  EDGE 76 77 -1  "
				"EDGE 143 134 -1  EDGE 143 152 -1" );
		reusable_FT_istream >> anchored_design_fold_tree;

		core::kinematics::FoldTree anchored_design_fold_tree_10at80;
		reusable_FT_istream.clear();
		reusable_FT_istream.str(
				"FOLD_TREE EDGE 1 59 -1  JEDGE 59 95 1 C N  INTRA_RES_STUB  EDGE 95 78 -1  EDGE 95 100 -1  "
				"EDGE 100 142 -1  EDGE 142 143 -1  EDGE 100 76 2  EDGE 142 153 3  EDGE 76 60 -1  EDGE 76 77 -1  "
				"EDGE 153 144 -1  EDGE 153 162 -1" );
		reusable_FT_istream >> anchored_design_fold_tree_10at80;

		TS_ASSERT_EQUALS( anchored_design_fold_tree_10at80,
				core::kinematics::remodel_fold_tree_to_account_for_insertion( anchored_design_fold_tree, 80, 10 ) );


		//Test complex case #2 (stolen from UBQ_E2_thioester_extra_bodies integration test)
		core::kinematics::FoldTree UBQ_E2_fold_tree;
		reusable_FT_istream.clear();
		reusable_FT_istream.str( "FOLD_TREE  EDGE 1 173 -1  EDGE 85 249 -2  SG   C    EDGE 249 174 -1  "
				"JEDGE 249 297 1 C NZ  END  EDGE 297 325 -1  EDGE 297 250 -1  JEDGE 1 326 2  N    N    END  "
				"EDGE 326 691 -1  JEDGE 1 692 3  N    N    END  EDGE 692 779 -1" );
		reusable_FT_istream >> UBQ_E2_fold_tree;

		core::kinematics::FoldTree UBQ_E2_fold_tree_10at300;
		reusable_FT_istream.clear();
		reusable_FT_istream.str(
				"FOLD_TREE  EDGE 1 173 -1  EDGE 85 249 -2  SG   C    EDGE 249 174 -1  JEDGE 249 297 1 C NZ  END  "
				"EDGE 297 335 -1  EDGE 297 250 -1  JEDGE 1 336 2  N    N    END  EDGE 336 701 -1  "
				"JEDGE 1 702 3  N    N    END  EDGE 702 789 -1" );
		reusable_FT_istream >> UBQ_E2_fold_tree_10at300;

		TS_ASSERT_EQUALS( UBQ_E2_fold_tree_10at300,
				core::kinematics::remodel_fold_tree_to_account_for_insertion( UBQ_E2_fold_tree, 300, 10 ) );
	}
};
