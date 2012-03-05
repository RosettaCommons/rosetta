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

	void test_more_complex_foldtree_setup(){

		core::pose::Pose pose;
		tr << "Test setup for multichain pose..."<< std::endl;

		//  DockingMultiChain.pdb
		//   A <- pdb numbering residues 1 through 214
		//   B <- pdb numbering residues 1 through 218
		//   E <- pdb numbering residues 1 through 129
		core::import_pose::pose_from_pdb(pose, "protocols/docking/DockingMultiChain.pdb" );
		core::kinematics::FoldTree fold_tree(pose.fold_tree());
		protocols::docking::DockJumps movable_jumps;

		//Specify partners via string:
		protocols::docking::setup_foldtree(pose, "AB_E", movable_jumps);
		std::stringstream partner_chainID_foldtree;
		partner_chainID_foldtree << pose.fold_tree();


		std::string database_filename("setup_docking_foldtree.db3");
		utility::file::file_delete(database_filename);
		utility::sql_database::sessionOP db_session(
			basic::database::get_db_session(
				database_filename, "sqlite3", false, true));

		std::string const schema(
			"CREATE TABLE IF NOT EXISTS interfaces (\n"
			"	struct_id INTEGER,\n"
			"	interface_id INTEGER,\n"
			"	PRIMARY KEY(struct_id, interface_id));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS interface_partners (\n"
			"	interface_id INTEGER,\n"
			"	partner_id INTEGER,\n"
			"	PRIMARY KEY(interface_id, partner_id));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS interface_partner_chains (\n"
			"	partner_id INTEGER,\n"
			"	chain_id INTEGER,\n"
			"	PRIMARY KEY(partner_id, chain_id));\n"
			"\n"
			"INSERT INTO interfaces VALUES (1, 1);\n"
			"INSERT INTO interface_partners VALUES (1, 1);\n"
			"INSERT INTO interface_partners VALUES (1, 2);\n"
			"INSERT INTO interface_partner_chains VALUES (1, 1);\n"
			"INSERT INTO interface_partner_chains VALUES (1, 2);\n"
			"INSERT INTO interface_partner_chains VALUES (2, 3);\n");

		basic::database::write_schema_to_database(schema, db_session);
		protocols::docking::setup_foldtree(pose, 1, db_session, movable_jumps);

		std::stringstream database_specified_foldtree;
		database_specified_foldtree << pose.fold_tree();

		TS_ASSERT_EQUALS(
			partner_chainID_foldtree.str(), database_specified_foldtree.str());

	}

};


