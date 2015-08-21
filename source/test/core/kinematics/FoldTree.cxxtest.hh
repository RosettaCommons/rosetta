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

//Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <sstream>

using core::kinematics::FoldTree;
using core::kinematics::Edge;

static thread_local basic::Tracer TR( "core.kinematics.FoldTree.cxxtest.hh" );

class FoldTreeTest : public CxxTest::TestSuite {
public:
	void setUp()
	{
		core_init();

		empty_ft_ = FoldTree();
		single_residue_ft_ = FoldTree( 1 );
		simple_ft_ = FoldTree( 40 );

		std::istringstream two_chain_ft_is, ligand_ft_is, single_loop_ft_is, reverse_ft_is, branching_ft_is,
			branching_ligand_ft_is;
		two_chain_ft_is.str( "FOLD_TREE  EDGE 1 20 -1  EDGE 1 21 1  EDGE 21 40 -1" );
		ligand_ft_is.str( "FOLD_TREE  EDGE 1 40 -1  EDGE 1 41 1" );
		single_loop_ft_is.str( "FOLD_TREE  EDGE 15 1 -1  EDGE 15 20 -1  EDGE 15 25 1  EDGE 25 21 -1  EDGE 25 40 -1" );
		reverse_ft_is.str( "FOLD_TREE  EDGE 40 1 -1" );
		branching_ft_is.str( "FOLD_TREE  EDGE 1 20 -1  EDGE 10 21 -2 X Z  EDGE 21 40 -1" );
		branching_ligand_ft_is.str( "FOLD_TREE  EDGE 1 40 -1  EDGE 1 41 1  EDGE 41 42 -1  EDGE 41 43 -2 X Z" );

		two_chain_ft_is >> two_chain_ft_;
		ligand_ft_is >> ligand_ft_;
		single_loop_ft_is >> single_loop_ft_;
		reverse_ft_is >> reverse_ft_;
		branching_ft_is >> branching_ft_;
		branching_ligand_ft_is >> branching_ligand_ft_;
	}

	void tearDown()
	{}

	void test_fold_tree_construction()
	{
		TS_ASSERT( ! empty_ft_.check_fold_tree() );
		TS_ASSERT( single_residue_ft_.check_fold_tree() );
		TS_ASSERT( simple_ft_.check_fold_tree() );
		TS_ASSERT( two_chain_ft_.check_fold_tree() );
		TS_ASSERT( ligand_ft_.check_fold_tree() );
		TS_ASSERT( single_loop_ft_.check_fold_tree() );
		TS_ASSERT( reverse_ft_.check_fold_tree() );
		TS_ASSERT( branching_ft_.check_fold_tree() );
		TS_ASSERT( branching_ligand_ft_.check_fold_tree() );

		std::ostringstream output;
		output << empty_ft_;
		TS_ASSERT_EQUALS( output.str(), "FOLD_TREE " );

		output.str( "" );
		output << single_residue_ft_;
		TS_ASSERT_EQUALS( output.str(), "FOLD_TREE  EDGE 1 1 -1 " );

		output.str( "" );
		output << simple_ft_;
		TS_ASSERT_EQUALS( output.str(), "FOLD_TREE  EDGE 1 40 -1 " );

		output.str( "" );
		output << two_chain_ft_;
		TS_ASSERT_EQUALS( output.str(), "FOLD_TREE  EDGE 1 20 -1  EDGE 1 21 1  EDGE 21 40 -1 " );

		output.str( "" );
		output << ligand_ft_;
		TS_ASSERT_EQUALS( output.str(), "FOLD_TREE  EDGE 1 40 -1  EDGE 1 41 1 " );

		output.str( "" );
		output << single_loop_ft_;
		TS_ASSERT_EQUALS( output.str(),
			"FOLD_TREE  EDGE 15 1 -1  EDGE 15 20 -1  EDGE 15 25 1  EDGE 25 21 -1  EDGE 25 40 -1 " );

		output.str( "" );
		output << reverse_ft_;
		TS_ASSERT_EQUALS( output.str(), "FOLD_TREE  EDGE 40 1 -1 " );

		output.str( "" );
		output << branching_ft_;
		TS_ASSERT_EQUALS( output.str(), "FOLD_TREE  EDGE 1 20 -1  EDGE 10 21 -2 X Z  EDGE 21 40 -1 " );

		output.str( "" );
		output << branching_ligand_ft_;
		TS_ASSERT_EQUALS( output.str(), "FOLD_TREE  EDGE 1 40 -1  EDGE 1 41 1  EDGE 41 42 -1  EDGE 41 43 -2 X Z " );
	}

	void test_boundary_right()
	{
		//TS_ASSERT_EQUALS( single_residue_ft_.boundary_right( 1 ), 1 );  // Why should this fail?

		// The latter function is a shortcut to the former.
		TS_ASSERT_EQUALS( simple_ft_.get_residue_edge( 21 ).stop(), simple_ft_.boundary_right( 21 ) );

		//TS_ASSERT_EQUALS( simple_ft_.boundary_right( 1 ), 40 );  // Why should this fail?
		TS_ASSERT_EQUALS( simple_ft_.boundary_right( 20 ), 40 );
		TS_ASSERT_EQUALS( simple_ft_.boundary_right( 40 ), 40 );

		TS_ASSERT_EQUALS( single_loop_ft_.boundary_right( 1 ), 1 );
		TS_ASSERT_EQUALS( single_loop_ft_.boundary_right( 10 ), 1 );
		//TS_ASSERT_EQUALS( single_loop_ft_.boundary_right( 15 ), 1 );  // Why should this fail?
		TS_ASSERT_EQUALS( single_loop_ft_.boundary_right( 20 ), 20 );
		TS_ASSERT_EQUALS( single_loop_ft_.boundary_right( 21 ), 21 );
		TS_ASSERT_EQUALS( single_loop_ft_.boundary_right( 25 ), 25 );
		TS_ASSERT_THROWS_ANYTHING( single_loop_ft_.boundary_right( 41 ) );  // not in tree!
	}

	void test_boundary_left()
	{
		//TS_ASSERT_EQUALS( single_residue_ft_.boundary_left( 1 ), 1 );  // Why should this fail?

		// The latter function is a shortcut to the former.
		TS_ASSERT_EQUALS( simple_ft_.get_residue_edge( 21 ).start(), simple_ft_.boundary_left( 21 ) );

		//TS_ASSERT_EQUALS( simple_ft_.boundary_left( 1 ), 1 );  // Why should this fail?
		TS_ASSERT_EQUALS( simple_ft_.boundary_left( 20 ), 1 );
		TS_ASSERT_EQUALS( simple_ft_.boundary_left( 40 ), 1 );

		TS_ASSERT_EQUALS( single_loop_ft_.boundary_left( 1 ), 15 );
		TS_ASSERT_EQUALS( single_loop_ft_.boundary_left( 10 ), 15 );
		//TS_ASSERT_EQUALS( single_loop_ft_.boundary_left( 15 ), 15 );  // Why should this fail?
		TS_ASSERT_EQUALS( single_loop_ft_.boundary_left( 20 ), 15 );
		TS_ASSERT_EQUALS( single_loop_ft_.boundary_left( 21 ), 25 );
		TS_ASSERT_EQUALS( single_loop_ft_.boundary_left( 25 ), 15 );
		TS_ASSERT_THROWS_ANYTHING( single_loop_ft_.boundary_left( 41 ) );  // not in tree!
	}


	// Ensure that repeated calls to hash value yield identical results when applied to the same FoldTree.
	void test_hash_value_unmodified() {
		TS_ASSERT_EQUALS( simple_ft_.hash_value(), simple_ft_.hash_value() );
	}

	// Ensure that operations that modify the FoldTree have an effect on the calculation of the hash value.
	void test_hash_value_modified() {
		FoldTree mod_ft( simple_ft_ );
		mod_ft.new_jump( 1, 3, 2 );

		TS_ASSERT_DIFFERS( simple_ft_.hash_value(), mod_ft.hash_value() );
	}

	// A quick check to see whether the FoldTree considers the final residues in a simple fold tree to be a cutpoint.
	void test_last_residue_is_cutpoint() {
		TS_ASSERT( simple_ft_.is_cutpoint( 40 ) );
	}

	void test_delete_self_edges()
	{
		TS_ASSERT( true );
	}

	void test_delete_seqpos()
	{
		TS_ASSERT( true );
	}

	void test_delete_jump_and_intervening_cutpoint()
	{
		TS_ASSERT( true );
	}

	void test_slide_cutpoint()
	{
		TS_ASSERT( true );
	}

	void test_slide_jump()
	{
		// create a simple fold tree with a jump
		std::stringstream ft_stream( "FOLD_TREE EDGE 1 2 1 EDGE 2 93 -1 " );
		core::kinematics::FoldTree ft;
		ft_stream >> ft;

		ft.slide_jump( 1, 1, 47 );
		TS_ASSERT( ft.check_fold_tree() );
		TS_ASSERT_EQUALS( ft.root(), 1 );

		core::kinematics::Edge e = ft.jump_edge( 1 );
		TS_ASSERT_EQUALS( (core::Size)e.start(), 1 );
		TS_ASSERT_EQUALS( (core::Size)e.stop(), 47 );
		TS_ASSERT_EQUALS( (core::Size)e.label(), 1 );
	}

	void test_delete_jump_seqpos()
	{
		TS_ASSERT( true );
	}

	void test_get_jump_that_builds_residue()
	{
		TS_ASSERT( true );
	}

	void test_get_parent_residue()
	{
		TS_ASSERT( true );
	}

	void test_insert_polymer_residue()
	{
		TS_ASSERT( true );
	}

	void test_insert_residue_by_chemical_bond()
	{
		TS_ASSERT( true );
	}

	void test_insert_residue_by_jump()
	{
		TS_ASSERT( true );
	}

	void test_insert_fold_tree_by_jump()
	{
		TS_ASSERT( true );
	}


	void test_append_residue()
	{
		TS_ASSERT( true );
	}

	void test_append_residue_by_chemical_bond()
	{
		TS_ASSERT( true );
	}

	void test_add_prepend_and_delete_edge_methods()
	{
		TS_ASSERT( true );
	}

	void test_delete_segment()
	{
		TS_ASSERT( true );
	}

	void test_update_edge_label()
	{
		TS_ASSERT( true );
	}

	void test_edge_label()
	{
		TS_ASSERT( true );
	}

	void test_delete_extra_vertices()
	{
		TS_ASSERT( true );
	}

	void test_downstream_and_upstream_jump_residue_methods()
	{
		TS_ASSERT( true );
	}

	void test_reorder()
	{
		TS_ASSERT( true );
	}

	void test_add_vertex()
	{
		TS_ASSERT( true );
	}

	void test_new_jump()
	{
		TS_ASSERT( true );
	}

	void test_new_chemical_bond()
	{
		TS_ASSERT( true );
	}

	void test_cut_edge()
	{
		TS_ASSERT( true );
	}

	void test_cutpoints()
	{
		TS_ASSERT( true );
	}

	void test_renumber_jumps()
	{
		TS_ASSERT( true );
	}

	void test_connected()
	{
		TS_ASSERT( true );
	}

	void test_partition_by_jump()
	{
		TS_ASSERT( true );
	}

	void test_partition_by_residue()
	{
		TS_ASSERT( true );
	}

	void test_cutpoint_by_jump()
	{
		TS_ASSERT( true );
	}

	void test_get_residue_edge()
	{
		TS_ASSERT( true );
	}

	void test_get_outgoing_edges()
	{
		TS_ASSERT( true );
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

	void test_get_polymer_residue_direction()
	{
		TS_ASSERT( true );
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

	void test_delete_jump() {
		core::kinematics::FoldTree ft;
		ft.add_edge( 37, 1, core::kinematics::Edge::PEPTIDE );
		ft.add_edge( 37, 74, core::kinematics::Edge::PEPTIDE );
		ft.add_edge( 37, 84, 1 );
		ft.add_edge( 84, 75, core::kinematics::Edge::PEPTIDE );
		ft.add_edge( 84, 93, 2 );
		ft.add_edge( 84, 92, core::kinematics::Edge::PEPTIDE );
		TS_ASSERT( ft.check_fold_tree() );
		ft.delete_jump_and_intervening_cutpoint( 1 );
		TS_ASSERT( ft.check_fold_tree() );
	}

private:
	FoldTree empty_ft_, single_residue_ft_, simple_ft_, two_chain_ft_, ligand_ft_, single_loop_ft_, reverse_ft_,
		branching_ft_, branching_ligand_ft_;
};
