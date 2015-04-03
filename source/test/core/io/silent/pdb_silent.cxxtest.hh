// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/silent/protein_silent.cxxtest.hh
/// @brief  test suite for protein silent-file format
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>


#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/import_pose/PDBSilentStruct.hh>
#include <core/scoring/rms_util.hh>

#include <utility/file/file_sys_util.hh>

//Auto Headers
#include <core/io/silent/EnergyNames.fwd.hh>
#include <utility/vector1.hh>


static basic::Tracer TR("test.core.io.silent.PDBSilentStruct");

using namespace core;

class PDBSilentStruct_Tests : public CxxTest::TestSuite
{
public:
	PDBSilentStruct_Tests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-in::file::silent_struct_type pdb" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}


void test_save_and_restore()
{
	double rms_threshold = 1e-3;
	pose::Pose start_pose(create_test_in_pdb_pose()), restored_pose;
	//core::import_pose::pose_from_pdb( start_pose, "core/io/test_in.pdb" );

	// Serialize the modified structure as a silent file
	core::chemical::ResidueTypeSetCOP	rsd_set
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	core::io::silent::SilentFileData sfd;
	std::string silent_outfile = "test_save_and_restore.test.silent.out";
	utility::file::file_delete( silent_outfile );
	core::import_pose::PDBSilentStruct pss( start_pose, "tag" );
	sfd.write_silent_struct( pss, silent_outfile );

	// Read the PDBSilentStruct from the silent-file
	sfd.read_file( silent_outfile );
	core::io::silent::SilentFileData::iterator iter = sfd.begin();
	TS_ASSERT( iter->decoy_tag() == "tag" );

	iter->fill_pose( restored_pose, *rsd_set );

	Real rms_to_restored = scoring::CA_rmsd( start_pose, restored_pose );
	TR << "RMS error from save/restore: " << rms_to_restored << std::endl;
	TS_ASSERT( rms_to_restored < rms_threshold );
	utility::file::file_delete( silent_outfile );


}

};
