// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/silent/binary_protein_silent_multipose.cxxtest.hh
/// @brief  Test suite for multiple poses stored in the same binary silent file.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/UTracer.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/VariantType.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <utility/file/file_sys_util.hh>

//Auto Headers
#include <core/io/silent/EnergyNames.fwd.hh>
#include <utility/vector1.hh>


static basic::Tracer TR("test.core.io.silent.protein_silent");

using namespace core;

class BinarySilentMultiposeTests : public CxxTest::TestSuite {

public:
	BinarySilentMultiposeTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-mute core.io.pdb -mute core.conformation -in::file::fullatom" );

		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCOP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		//ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		//  if(!residue_set.has_name("GTP")) params_files.push_back("core/io/GTP.params");
		//  residue_set.read_files(params_files);
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_save_multiple_poses_and_restore()
	{
		TR << "Starting test_save_multiple_poses_and_restore()." << std::endl;
		TR << "Author: Vikram K. Mulligan, Baker lab (vmullig@uw.edu)." << std::endl;
		TR << "If this test were to fail, it would indicate that binary silent files are not properly storing multiple poses." << std::endl;

		double rms_threshold = 1e-2;
		double score_threshold = 1e-1;

		pose::Pose ref_pose, restored_pose, ref_pose_2, restored_pose_2;
		core::chemical::ResidueTypeSetCOP rsd =
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

		//Import the poses:
		core::import_pose::pose_from_pdb( ref_pose, *rsd, std::string("core/io/test_in.pdb"));
		core::import_pose::pose_from_pdb( ref_pose_2, *rsd, std::string("core/io/bin_silentfile_test.pdb"));

		//Write out the poses:
		std::string const silent_outfile( "core/io/bin_silentfile_multipose_test.out" );
		utility::file::file_delete( silent_outfile );
		core::io::silent::SilentFileData sfd;
		core::io::silent::BinarySilentStruct pss( ref_pose, "tag" );
		sfd.write_silent_struct( pss, silent_outfile );
		core::io::silent::BinarySilentStruct pss2( ref_pose_2, "tag2" );
		sfd.write_silent_struct( pss2, silent_outfile );


		// Read the ProteinSilentStruct from the silent-file
		core::io::silent::SilentFileData sfd2;
		sfd2.read_file( silent_outfile );
		TS_ASSERT( sfd2.size() > 0 );
		core::io::silent::SilentFileData::iterator iter = sfd2.begin();
		iter->fill_pose( restored_pose, *rsd );
		++iter;
		iter->fill_pose( restored_pose_2, *rsd );

		// test rms difference
		Real rms_to_restored = scoring::CA_rmsd( ref_pose, restored_pose );
		TR << "RMS error from save/restore #1: " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < rms_threshold );

		Real rms_to_restored_2 = scoring::CA_rmsd( ref_pose_2, restored_pose_2 );
		TR << "RMS error from save/restore #2: " << rms_to_restored_2 << std::endl;
		TS_ASSERT( rms_to_restored_2 < rms_threshold );


		// test score13 difference
		core::scoring::ScoreFunctionOP scorefxn =
			core::scoring::ScoreFunctionFactory::create_score_function( "score13_env_hb" );
		Real score_ref = (*scorefxn)(ref_pose);
		Real score_restored = (*scorefxn)(restored_pose);
		Real score_del = std::fabs( score_restored - score_ref );
		TR << "Score difference #1: " << score_del << std::endl;
		if ( score_del > score_threshold ) {
			restored_pose.dump_pdb( "restored_pose_wrong_score.pdb" );
			ref_pose.dump_pdb( "ref_pose_wrong_score.pdb" );
		}
		TS_ASSERT( score_del < score_threshold );

		Real score_ref_2 = (*scorefxn)(ref_pose_2);
		Real score_restored_2 = (*scorefxn)(restored_pose_2);
		Real score_del_2 = std::fabs( score_restored_2 - score_ref_2 );
		TR << "Score difference #2: " << score_del_2 << std::endl;
		if ( score_del_2 > score_threshold ) {
			restored_pose_2.dump_pdb( "restored_pose_wrong_score_2.pdb" );
			ref_pose_2.dump_pdb( "ref_pose_wrong_score_2.pdb" );
		}
		TS_ASSERT( score_del_2 < score_threshold );

		// utility::file::file_delete( silent_outfile );

	}

	void test_save_multiple_poses_with_resize_and_restore()
	{
		TR << "Starting test_save_multiple_poses_with_resize_and_restore()." << std::endl;
		TR << "Author: Vikram K. Mulligan, Baker lab (vmullig@uw.edu)." << std::endl;
		TR << "If this test were to fail, it would indicate that binary silent files are not properly storing multiple poses when one of the poses is resized." << std::endl;

		double rms_threshold = 1e-2;
		double score_threshold = 1e-1;

		pose::Pose ref_pose, restored_pose, ref_pose_2, restored_pose_2;
		core::chemical::ResidueTypeSetCOP rsd =
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		core::import_pose::pose_from_pdb( ref_pose, *rsd, std::string("core/io/test_in.pdb"));
		core::import_pose::pose_from_pdb( ref_pose_2, *rsd, std::string("core/io/bin_silentfile_test.pdb"));

		//Add a residue to ref_pose_2 ONLY (and leave ref_pose alone!)
		core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( rsd->name_map( "ALA" ) ) );
		core::Size nres = ref_pose_2.n_residue();
		if ( ref_pose_2.residue(nres).has_variant_type(core::chemical::UPPER_TERMINUS_VARIANT) ) {
			core::pose::remove_variant_type_from_pose_residue( ref_pose_2, core::chemical::UPPER_TERMINUS_VARIANT, nres );
		}
		ref_pose_2.append_residue_by_bond(*new_rsd, true, 0, nres, 0, false);
		core::pose::add_variant_type_to_pose_residue( ref_pose_2, core::chemical::UPPER_TERMINUS_VARIANT, ref_pose_2.n_residue() );

		std::string const silent_outfile( "core/io/bin_silentfile_multipose_test.out" );
		utility::file::file_delete( silent_outfile );
		core::io::silent::SilentFileData sfd;
		core::io::silent::BinarySilentStruct pss( ref_pose, "tag" );
		sfd.write_silent_struct( pss, silent_outfile );
		core::io::silent::BinarySilentStruct pss2( ref_pose_2, "tag2" );
		sfd.write_silent_struct( pss2, silent_outfile );


		// Read the ProteinSilentStruct from the silent-file
		core::io::silent::SilentFileData sfd2;
		sfd2.read_file( silent_outfile );
		TS_ASSERT( sfd2.size() > 0 );
		core::io::silent::SilentFileData::iterator iter = sfd2.begin();
		iter->fill_pose( restored_pose, *rsd );
		++iter;
		iter->fill_pose( restored_pose_2, *rsd );

		// test rms difference
		Real rms_to_restored = scoring::CA_rmsd( ref_pose, restored_pose );
		TR << "RMS error from save/restore #1: " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < rms_threshold );

		Real rms_to_restored_2 = scoring::CA_rmsd( ref_pose_2, restored_pose_2 );
		TR << "RMS error from save/restore #2: " << rms_to_restored_2 << std::endl;
		TS_ASSERT( rms_to_restored_2 < rms_threshold );


		// test score13 difference
		core::scoring::ScoreFunctionOP scorefxn =
			core::scoring::ScoreFunctionFactory::create_score_function( "score13_env_hb" );
		Real score_ref = (*scorefxn)(ref_pose);
		Real score_restored = (*scorefxn)(restored_pose);
		Real score_del = std::fabs( score_restored - score_ref );
		TR << "Score difference #1: " << score_del << std::endl;
		if ( score_del > score_threshold ) {
			restored_pose.dump_pdb( "restored_pose_wrong_score.pdb" );
			ref_pose.dump_pdb( "ref_pose_wrong_score.pdb" );
		}
		TS_ASSERT( score_del < score_threshold );

		Real score_ref_2 = (*scorefxn)(ref_pose_2);
		Real score_restored_2 = (*scorefxn)(restored_pose_2);
		Real score_del_2 = std::fabs( score_restored_2 - score_ref_2 );
		TR << "Score difference #2: " << score_del_2 << std::endl;
		if ( score_del_2 > score_threshold ) {
			restored_pose_2.dump_pdb( "restored_pose_wrong_score_2.pdb" );
			ref_pose_2.dump_pdb( "ref_pose_wrong_score_2.pdb" );
		}
		TS_ASSERT( score_del_2 < score_threshold );

		// utility::file::file_delete( silent_outfile );

	}

};
