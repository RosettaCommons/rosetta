// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/PDB_w_lipids_IO.cxxtest.hh
/// @brief  test suite for reader/writer for PDBs with lipids
/// @author Rocco Moretti

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <core/io/pose_from_sfr/PoseFromSFRBuilder.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/FoldTree.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("core.io.PDB_w_lipids_IO.cxxtest");

using namespace core;

class lipid_PDBs_IO : public CxxTest::TestSuite
{

public:
	lipid_PDBs_IO() {}

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-include_lipids" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_branched() {
		pose::Pose pose;
		std::string const original_file_name("core/io/two_lipids.pdb");
		import_pose::pose_from_file(pose, original_file_name, core::import_pose::PDB_file);

		TS_ASSERT_EQUALS( pose.size(), 8 );
		TS_ASSERT_EQUALS( pose.residue_type(1).name(), "glycerol-3-phosphate" );
		TS_ASSERT_EQUALS( pose.residue_type(2).name(), "choline" );
		TS_ASSERT_EQUALS( pose.residue_type(3).name(), "palmitate" );
		TS_ASSERT_EQUALS( pose.residue_type(4).name(), "oleate" );
		TS_ASSERT_EQUALS( pose.residue_type(5).name(), "glycerol-3-phosphate" );
		TS_ASSERT_EQUALS( pose.residue_type(6).name(), "choline" );
		TS_ASSERT_EQUALS( pose.residue_type(7).name(), "palmitate" );
		TS_ASSERT_EQUALS( pose.residue_type(8).name(), "oleate" );

		TS_ASSERT( pose.residue(1).is_lower_terminus() );
		TS_ASSERT( pose.residue(1).is_branch_point() );
		TS_ASSERT( pose.residue(2).is_upper_terminus() );
		TS_ASSERT( pose.residue(3).is_upper_terminus() );
		TS_ASSERT( pose.residue(4).is_upper_terminus() );
		TS_ASSERT( pose.residue(5).is_lower_terminus() );
		TS_ASSERT( pose.residue(5).is_branch_point() );
		TS_ASSERT( pose.residue(6).is_upper_terminus() );
		TS_ASSERT( pose.residue(7).is_upper_terminus() );
		TS_ASSERT( pose.residue(8).is_upper_terminus() );

		TS_ASSERT_EQUALS( pose.residue(2).connected_residue_at_lower(), 1 );
		TS_ASSERT_EQUALS( pose.residue(3).connected_residue_at_lower(), 1 );
		TS_ASSERT_EQUALS( pose.residue(4).connected_residue_at_lower(), 1 );
		TS_ASSERT_EQUALS( pose.residue(6).connected_residue_at_lower(), 5 );
		TS_ASSERT_EQUALS( pose.residue(7).connected_residue_at_lower(), 5 );
		TS_ASSERT_EQUALS( pose.residue(8).connected_residue_at_lower(), 5 );

		TS_ASSERT_EQUALS( pose.residue(1).connected_residue_at_upper(), 2 );
		TS_ASSERT_EQUALS( pose.residue(5).connected_residue_at_upper(), 6 );

		TS_ASSERT_EQUALS( pose.residue(1).residue_connection_partner(2), 3 );
		TS_ASSERT_EQUALS( pose.residue(1).residue_connection_partner(3), 4 );
		TS_ASSERT_EQUALS( pose.residue(5).residue_connection_partner(2), 7 );
		TS_ASSERT_EQUALS( pose.residue(5).residue_connection_partner(3), 8 );

		TS_ASSERT_EQUALS( pose.fold_tree().to_string(), "FOLD_TREE  EDGE 1 2 -1  EDGE 1 3 -2  O1   C1   EDGE 1 4 -2  O2   C1   EDGE 1 5 1  EDGE 5 6 -1  EDGE 5 7 -2  O1   C1   EDGE 5 8 -2  O2   C1  " )
	}

};

