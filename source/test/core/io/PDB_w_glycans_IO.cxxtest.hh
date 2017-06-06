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

// Project Headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("core.io.PDB_w_glycans_IO.cxxtest");

using namespace core;

class PDB_w_glycans_IO : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCAP residue_set;

public:
	PDB_w_glycans_IO() {}

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-auto_detect_glycan_connections -include_sugars -alternate_3_letter_codes pdb_sugar -ignore_zero_occupancy -s 5FYL.pdb -write_pdb_link_records" );
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// ------------------------------------------ //
	// test that PDB input and output function properly
	void test_pdb_io() {
		pose::Pose pose;
		const std::string original_file_name("core/io/5FYL.pdb");
		import_pose::pose_from_file(pose, original_file_name, core::import_pose::PDB_file);

		// write pose to a new file...
		const std::string tmp_file_name("PDB_w_glycans_IO_cxxtest.pdb._tmp_");
		io::pdb::dump_pdb(pose, tmp_file_name);

		const std::string correct_output_file_name("core/io/5FYL_correct_output.pdb");
		TS_ASSERT_FILE_EQ(correct_output_file_name.c_str(), tmp_file_name.c_str());

	}

};
