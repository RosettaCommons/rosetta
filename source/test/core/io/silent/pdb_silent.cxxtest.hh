// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/silent/protein_silent.cxxtest.hh
/// @brief  test suite for protein silent-file format
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>


#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/import_pose/PDBSilentStruct.hh>
#include <core/scoring/rms_util.hh>
#include <core/chemical/ResidueTypeSet.hh>

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
		core_init_with_additional_options( "-in::file::silent_struct_type pdb -out:level 500" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_save_and_restore()
	{
		double rms_threshold = 1e-3;
		pose::Pose start_pose(create_test_in_pdb_pose()), restored_pose;
		//core::import_pose::pose_from_file( start_pose, "core/io/test_in.pdb" , core::import_pose::PDB_file);

		// Serialize the modified structure as a silent file
		core::chemical::ResidueTypeSetCOP rsd_set
			= core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

		core::io::silent::SilentFileOptions opts;
		core::io::silent::SilentFileData sfd( opts );
		std::string silent_outfile = "test_save_and_restore.test.silent.out";
		utility::file::file_delete( silent_outfile );
		core::import_pose::PDBSilentStruct pss( opts, start_pose, "tag" );
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

	void test_pdb_silent_disulfide() {
		pose::Pose start_pose(fullatom_pose_from_string(disulfide_string())), restored_pose;
		core::chemical::ResidueTypeSetCOP rsd_set
			= core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

		core::chemical::ResidueType const cydn = rsd_set->name_map( "CYS:NtermProteinFull:disulfide" );
		core::chemical::ResidueType const cydc = rsd_set->name_map( "CYS:CtermProteinFull:disulfide" );
		TS_ASSERT( start_pose.residue( 1 ).type().name() == cydn.name() );
		TS_ASSERT( start_pose.residue( 2 ).type().name() == cydc.name() );

		core::io::silent::SilentFileOptions opts;
		core::io::silent::SilentFileData sfd( opts );
		std::string silent_outfile = "test.out";
		utility::file::file_delete( silent_outfile );
		core::import_pose::PDBSilentStruct pss( opts, start_pose, "test_0001" );
		sfd.write_silent_struct( pss, silent_outfile );


		// Read the PDBSilentStruct from the silent-file
		sfd.read_file( silent_outfile );
		core::io::silent::SilentFileData::iterator iter = sfd.begin();
		TS_ASSERT( iter->decoy_tag() == "test_0001" );

		iter->fill_pose( restored_pose, *rsd_set );

		// can't test type equality, can test names?
		TS_ASSERT( restored_pose.residue( 1 ).type().name() == cydn.name() );
		TS_ASSERT( restored_pose.residue( 2 ).type().name() == cydc.name() );

		// round trip successful
		// now test compatibility with how Old Rosetta might render it
		// This one reads in two residues that have BOTH termini
		// the reason it reads in differently is because earlier, we had done pose from string
		// and then written out before any actions that would have completed termini

		core::chemical::ResidueType const cyd = rsd_set->name_map( "CYS:CtermProteinFull:NtermProteinFull:disulfide" );
		core::io::silent::SilentFileData sfd2( opts );
		std::string silent_outfile2 = "core/io/test_CYD.out";

		// Read the PDBSilentStruct from the silent-file
		sfd2.read_file( silent_outfile2 );
		core::io::silent::SilentFileData::iterator iter2 = sfd2.begin();
		TS_ASSERT( iter2->decoy_tag() == "test_0001" );

		iter2->fill_pose( restored_pose, *rsd_set );

		TS_ASSERT( restored_pose.residue( 1 ).type().name() == cyd.name() );
		TS_ASSERT( restored_pose.residue( 2 ).type().name() == cyd.name() );

	}
	inline
	std::string
	disulfide_string()
	{
		return
			"ATOM      1  N   CYS F   6      67.701 112.472  49.297  1.00  0.00           N\n"
			"ATOM      2  CA  CYS F   6      68.061 111.684  50.465  1.00  0.00           C\n"
			"ATOM      3  C   CYS F   6      67.597 110.245  50.279  1.00  0.00           C\n"
			"ATOM      4  O   CYS F   6      68.213 109.307  50.783  1.00  0.00           O\n"
			"ATOM      5  OXT CYS F   6      66.614 110.010  49.631  1.00  0.00           O\n"
			"ATOM      6  CB  CYS F   6      67.288 112.366  51.593  1.00  0.00           C\n"
			"ATOM      7  SG  CYS F   6      67.791 114.073  51.920  1.00  0.00           S\n"
			"ATOM      8 1H   CYS F   6      68.007 113.415  49.423  1.00  0.00           H\n"
			"ATOM      9 2H   CYS F   6      68.139 112.086  48.486  1.00  0.00           H\n"
			"ATOM     10 3H   CYS F   6      66.708 112.458  49.178  1.00  0.00           H\n"
			"ATOM     11  HA  CYS F   6      69.135 111.701  50.568  1.00  0.00           H\n"
			"ATOM     12 1HB  CYS F   6      66.224 112.405  51.352  1.00  0.00           H\n"
			"ATOM     13 2HB  CYS F   6      67.432 111.825  52.528  1.00  0.00           H\n"
			"ATOM     14  N   CYS F  19      72.082 115.921  51.930  1.00  0.00           N\n"
			"ATOM     15  CA  CYS F  19      71.424 115.188  52.997  1.00  0.00           C\n"
			"ATOM     16  C   CYS F  19      72.429 114.816  54.084  1.00  0.00           C\n"
			"ATOM     17  O   CYS F  19      73.633 114.718  53.816  1.00  0.00           O\n"
			"ATOM     18  OXT CYS F  19      72.052 114.614  55.205  1.00  0.00           O\n"
			"ATOM     19  CB  CYS F  19      70.907 113.938  52.285  1.00  0.00           C\n"
			"ATOM     20  SG  CYS F  19      69.671 114.264  51.005  1.00  0.00           S\n"
			"ATOM     21 1H   CYS F  19      71.413 116.158  51.226  1.00  0.00           H\n"
			"ATOM     22 2H   CYS F  19      72.486 116.757  52.300  1.00  0.00           H\n"
			"ATOM     23 3H   CYS F  19      72.801 115.353  51.529  1.00  0.00           H\n"
			"ATOM     24  HA  CYS F  19      70.664 115.822  53.422  1.00  0.00           H\n"
			"ATOM     25 1HB  CYS F  19      71.727 113.419  51.788  1.00  0.00           H\n"
			"ATOM     26 2HB  CYS F  19      70.431 113.266  53.000  1.00  0.00           H\n"
			"TER\n";
	}


	/*inline
	std::string
	homocysteine_disulfide_string()
	{
	return
	"ATOM      1  N   CYS F   6      67.701 112.472  49.297  1.00  0.00           N\n"
	"ATOM      2  CA  CYS F   6      68.061 111.684  50.465  1.00  0.00           C\n"
	"ATOM      3  C   CYS F   6      67.597 110.245  50.279  1.00  0.00           C\n"
	"ATOM      4  O   CYS F   6      68.213 109.307  50.783  1.00  0.00           O\n"
	"ATOM      5  OXT CYS F   6      66.614 110.010  49.631  1.00  0.00           O\n"
	"ATOM      6  CB  CYS F   6      67.288 112.366  51.593  1.00  0.00           C\n"
	"ATOM      7  SG  CYS F   6      67.791 114.073  51.920  1.00  0.00           S\n"
	"ATOM      8 1H   CYS F   6      68.007 113.415  49.423  1.00  0.00           H\n"
	"ATOM      9 2H   CYS F   6      68.139 112.086  48.486  1.00  0.00           H\n"
	"ATOM     10 3H   CYS F   6      66.708 112.458  49.178  1.00  0.00           H\n"
	"ATOM     11  HA  CYS F   6      69.135 111.701  50.568  1.00  0.00           H\n"
	"ATOM     12 1HB  CYS F   6      66.224 112.405  51.352  1.00  0.00           H\n"
	"ATOM     13 2HB  CYS F   6      67.432 111.825  52.528  1.00  0.00           H\n"
	"ATOM     14  N   CYS F  19      72.082 115.921  51.930  1.00  0.00           N\n"
	"ATOM     15  CA  CYS F  19      71.424 115.188  52.997  1.00  0.00           C\n"
	"ATOM     16  C   CYS F  19      72.429 114.816  54.084  1.00  0.00           C\n"
	"ATOM     17  O   CYS F  19      73.633 114.718  53.816  1.00  0.00           O\n"
	"ATOM     18  OXT CYS F  19      72.052 114.614  55.205  1.00  0.00           O\n"
	"ATOM     19  CB  CYS F  19      70.907 113.938  52.285  1.00  0.00           C\n"
	"ATOM     20  SG  CYS F  19      69.671 114.264  51.005  1.00  0.00           S\n"
	"ATOM     21 1H   CYS F  19      71.413 116.158  51.226  1.00  0.00           H\n"
	"ATOM     22 2H   CYS F  19      72.486 116.757  52.300  1.00  0.00           H\n"
	"ATOM     23 3H   CYS F  19      72.801 115.353  51.529  1.00  0.00           H\n"
	"ATOM     24  HA  CYS F  19      70.664 115.822  53.422  1.00  0.00           H\n"
	"ATOM     25 1HB  CYS F  19      71.727 113.419  51.788  1.00  0.00           H\n"
	"ATOM     26 2HB  CYS F  19      70.431 113.266  53.000  1.00  0.00           H\n"
	"TER\n";
	}*/
};


