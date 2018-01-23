// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/PDB_w_glycans_IO.cxxtest.hh
/// @brief  test suite for reader/writer for PDBs with glycans
/// @author Sebastian Rämisch

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

static basic::Tracer TR("core.io.PDB_w_glycans_IO.cxxtest");

using namespace core;

class glycan_PDBs_IO : public CxxTest::TestSuite
{

public:
	glycan_PDBs_IO() {}

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-include_sugars -output_virtual" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_branched() {
		pose::Pose pose;
		std::string const original_file_name("core/chemical/carbohydrates/amylopectin_fragment.pdb");
		import_pose::pose_from_file(pose, original_file_name, core::import_pose::PDB_file);

		TS_ASSERT_EQUALS( pose.size(), 6 );
		TS_ASSERT_EQUALS( pose.residue_type(1).name(), "->4)-alpha-D-Glcp:reducing_end" );
		TS_ASSERT_EQUALS( pose.residue_type(2).name(), "->4)-alpha-D-Glcp:->6)-branch" );
		TS_ASSERT_EQUALS( pose.residue_type(3).name(), "->4)-alpha-D-Glcp" );
		TS_ASSERT_EQUALS( pose.residue_type(4).name(), "->4)-alpha-D-Glcp:non-reducing_end" );
		TS_ASSERT_EQUALS( pose.residue_type(5).name(), "->4)-alpha-D-Glcp" );
		TS_ASSERT_EQUALS( pose.residue_type(6).name(), "->4)-alpha-D-Glcp:non-reducing_end" );

		TS_ASSERT( pose.residue(1).is_lower_terminus() );
		TS_ASSERT( pose.residue(2).is_branch_point() );
		TS_ASSERT( pose.residue(4).is_upper_terminus() );
		TS_ASSERT( pose.residue(6).is_upper_terminus() );

		TS_ASSERT_EQUALS( pose.residue(2).connected_residue_at_lower(), 1 );
		TS_ASSERT_EQUALS( pose.residue(3).connected_residue_at_lower(), 2 );
		TS_ASSERT_EQUALS( pose.residue(4).connected_residue_at_lower(), 3 );
		TS_ASSERT_EQUALS( pose.residue(5).connected_residue_at_lower(), 2 );
		TS_ASSERT_EQUALS( pose.residue(6).connected_residue_at_lower(), 5 );
		TS_ASSERT_EQUALS( pose.residue(1).connected_residue_at_upper(), 2 );
		TS_ASSERT_EQUALS( pose.residue(2).connected_residue_at_upper(), 3 );
		TS_ASSERT_EQUALS( pose.residue(3).connected_residue_at_upper(), 4 );
		TS_ASSERT_EQUALS( pose.residue(5).connected_residue_at_upper(), 6 );
		TS_ASSERT_EQUALS( pose.residue(2).residue_connection_partner(3), 5 );
		TS_ASSERT_EQUALS( pose.fold_tree().to_string(), "FOLD_TREE  EDGE 1 4 -1  EDGE 2 5 -2  O6   C1   EDGE 5 6 -1 " )
	}

};

// Separate class for `-auto_detect_glycan_connections`
class PDB_w_glycans_IO : public CxxTest::TestSuite
{

public:
	PDB_w_glycans_IO() {}

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-auto_detect_glycan_connections -include_sugars -alternate_3_letter_codes pdb_sugar -ignore_zero_occupancy -write_pdb_link_records" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_glycan_fragment() {
		pose::Pose pose;
		import_pose::pose_from_file(pose, "core/io/5FYL_frag_1.pdb", core::import_pose::PDB_file);

		TS_ASSERT_EQUALS( pose.size(), 55 );
		TS_ASSERT_EQUALS( pose.residue(11).name(), "ASN:N-glycosylated" );
		TS_ASSERT_EQUALS( pose.residue(18).name(), "ASN:N-glycosylated" );
		TS_ASSERT_EQUALS( pose.residue(37).name(), "ASN:N-glycosylated" );
		TS_ASSERT_EQUALS( pose.residue(50).name(), "GLN:CtermProteinFull" );
		TS_ASSERT_EQUALS( pose.residue(51).name(), "->3)-beta-D-Glcp:non-reducing_end:2-AcNH" );
		TS_ASSERT_EQUALS( pose.residue(52).name(), "->3)-beta-D-Glcp:non-reducing_end:2-AcNH" );
		TS_ASSERT_EQUALS( pose.residue(53).name(), "->4)-beta-D-Glcp:2-AcNH" );
		TS_ASSERT_EQUALS( pose.residue(54).name(), "->4)-beta-D-Glcp:2-AcNH" );
		TS_ASSERT_EQUALS( pose.residue(55).name(), "->3)-beta-D-Manp:non-reducing_end" );

		TS_ASSERT_EQUALS( pose.residue(51).connected_residue_at_lower(), 11 );
		TS_ASSERT_EQUALS( pose.residue(52).connected_residue_at_lower(), 18 );
		TS_ASSERT_EQUALS( pose.residue(53).connected_residue_at_lower(), 37 );
		TS_ASSERT_EQUALS( pose.residue(54).connected_residue_at_lower(), 53 );
		TS_ASSERT_EQUALS( pose.residue(55).connected_residue_at_lower(), 54 );

		TS_ASSERT_EQUALS( pose.residue(53).connected_residue_at_upper(), 54 );
		TS_ASSERT_EQUALS( pose.residue(54).connected_residue_at_upper(), 55 );

		TS_ASSERT_EQUALS( pose.fold_tree().to_string(), "FOLD_TREE  EDGE 1 50 -1  EDGE 11 51 -2  ND2  C1   EDGE 18 52 -2  ND2  C1   EDGE 37 53 -2  ND2  C1   EDGE 53 55 -1 " )
	}

	// ------------------------------------------ //
	// test that PDB input and output function properly
	void test_pdb_io() {
		pose::Pose pose;
		std::string const original_file_name("core/io/5FYL.pdb");
		import_pose::pose_from_file(pose, original_file_name, core::import_pose::PDB_file);

		// Consistent scoring (but with a score table in output.)
		core::scoring::ScoreFunction empty_sfxn;
		empty_sfxn(pose);

		// write pose to a new file...
		std::string const tmp_file_name("PDB_w_glycans_IO_cxxtest.pdb._tmp_");
		io::pdb::dump_pdb(pose, tmp_file_name);

		TR << "5FYL FT: " << pose.fold_tree() << std::endl;

		///Make sure we can load it back in after export.
		TS_ASSERT_THROWS_NOTHING(import_pose::pose_from_file(pose,tmp_file_name));

		//Test the contents.
		std::string const correct_output_file_name("core/io/5FYL_correct_output.pdb");
		TS_ASSERT_FILE_EQ(correct_output_file_name.c_str(), tmp_file_name.c_str());



	}

};
