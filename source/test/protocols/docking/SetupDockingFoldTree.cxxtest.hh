// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/docking/SetupDockingFoldTree.cxxtest.hh
/// @brief  tests for setting up fold trees across interfaces by specifing the chains of the partners
/// @author Matthew O'Meara mattjomeara@gmail.com

// Unit headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <protocols/docking/util.hh>


// project headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/database/sql_utils.hh>

#include <basic/Tracer.hh>

#include <sstream>

static basic::Tracer tr("protocols.docking.SetupDockingFoldTree.cxxtest");

class SetupDockingFoldTree : public CxxTest::TestSuite {

public:
	void setUp() {
		core_init();
	}

	void tearDown() {}

	/// @brief test the docking protocol functions
	void test_autodetct_vs_specifing_first_jump() {

		tr << "Test setup of simple foldtree..." << std::endl;

		core::pose::Pose pose;

		// DockingTest.pdb has two chains
		//  E <- pdb numbering residues 1 through 245
		//  I <- pdb numbering residues 1 through 56
		core::import_pose::pose_from_pdb(pose, "protocols/docking/DockingTest.pdb");
		core::Size const rb_jump(1);

		//setting up the fold tree as is used in docking
		core::kinematics::FoldTree fold_tree;

		core::Size jump_pos1 = 197;
		core::Size jump_pos2 = 282;
		core::Size cutpoint = 245;

		fold_tree.clear();
		fold_tree.add_edge( jump_pos1, jump_pos2, rb_jump );
		fold_tree.add_edge( 1, jump_pos1, core::kinematics::Edge::PEPTIDE );
		fold_tree.add_edge( jump_pos1, cutpoint, core::kinematics::Edge::PEPTIDE );
		fold_tree.add_edge( jump_pos2, cutpoint+1, core::kinematics::Edge::PEPTIDE );
		fold_tree.add_edge( jump_pos2, pose.total_residue(), core::kinematics::Edge::PEPTIDE );
		fold_tree.reorder( 1 );

		pose.fold_tree(fold_tree);
		protocols::docking::DockJumps movable_jumps;

		tr << "Initial foldtree:" << std::endl;
		tr << "\t" << pose.fold_tree() << std::endl;
		std::stringstream default_foldtree, specified_foldtree;

		//autodetect first jump:
		protocols::docking::setup_foldtree(pose, "_", movable_jumps);
		default_foldtree << pose.fold_tree();
		pose.fold_tree(fold_tree);

		//specify first jump via string:
		protocols::docking::setup_foldtree(pose, "E_I", movable_jumps);
		specified_foldtree << pose.fold_tree();
		pose.fold_tree(fold_tree);

		TS_ASSERT_EQUALS(default_foldtree.str(), specified_foldtree.str());
	}
	
	void test_residue_selector_based_setup(){
		core::pose::Pose pose;
		tr << "Test setup for multichain pose..."<< std::endl;
		core::import_pose::pose_from_pdb(pose, "protocols/docking/DockingMultiChain.pdb" );
		protocols::docking::DockJumps movable_jumps;
		movable_jumps.push_back( 1 );
		
		std::stringstream target_foldtree;
		target_foldtree << "FOLD_TREE  ";
		target_foldtree << "EDGE 1 214 -1  ";
		target_foldtree << "EDGE 214 215 2  ";
		target_foldtree << "EDGE 215 384 -1  ";
		target_foldtree << "EDGE 384 432 -1  ";
		target_foldtree << "EDGE 384 488 1  ";
		target_foldtree << "EDGE 488 433 -1  ";
		target_foldtree << "EDGE 488 561 -1 ";
		
		protocols::docking::setup_foldtree( pose, "AB_E", movable_jumps );
		
		std::stringstream computed_foldtree;
		computed_foldtree << pose.fold_tree();
		
		TS_ASSERT_EQUALS( computed_foldtree.str(), target_foldtree.str() );
		
		std::stringstream out_of_order_foldtree;
		out_of_order_foldtree << "FOLD_TREE  ";
		out_of_order_foldtree << "EDGE 1 98 -1  ";
		out_of_order_foldtree << "EDGE 98 214 -1  ";
		out_of_order_foldtree << "EDGE 98 325 1  ";
		out_of_order_foldtree << "EDGE 214 433 2  ";
		out_of_order_foldtree << "EDGE 325 215 -1  ";
		out_of_order_foldtree << "EDGE 325 432 -1  ";
		out_of_order_foldtree << "EDGE 433 561 -1 ";
		
		protocols::docking::setup_foldtree( pose, "AE_B", movable_jumps );

		std::stringstream computed_out_of_order_foldtree;
		computed_out_of_order_foldtree << pose.fold_tree();
		
		TS_ASSERT_EQUALS( computed_out_of_order_foldtree.str(), out_of_order_foldtree.str() );
	}

};


