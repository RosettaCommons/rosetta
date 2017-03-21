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
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- added test for badly-formatted silent files of a sort commonly produced by Rosetta.

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
#include <core/io/silent/SilentFileOptions.hh>
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

class BinarySilentTests : public CxxTest::TestSuite {

public:
	BinarySilentTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-mute core.io.pdb -mute core.conformation -in::file::fullatom" );
		//core_init_with_additional_options( "-out:level 4000 -in::file::fullatom" );

		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCOP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		//ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		//  if(!residue_set.has_name("GTP")) params_files.push_back("core/io/GTP.params");
		//  residue_set.read_files_for_custom_residue_types(params_files);
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_save_and_restore_centroid() {
		using namespace core::chemical;
		ResidueTypeSetCOP cen_rsd_set =
			core::chemical::ChemicalManager::get_instance()->residue_type_set(
			"centroid"
		);
		pose::Pose ref_pose, restored_pose;
		core::import_pose::pose_from_file( ref_pose, *cen_rsd_set, "core/io/bin_silentfile_test.pdb", core::import_pose::PDB_file);
		double rms_threshold = 1e-2;
		double score_threshold = 1e-1;
		TS_ASSERT( !ref_pose.is_fullatom() );
		// Read the ProteinSilentStruct from the silent-file
		core::io::silent::SilentFileOptions opts;
		core::io::silent::SilentFileData sfd( opts );
		std::string const silent_outfile( "core/io/bin_silentfile_centroid.out" ); // read file w/ non-ideal geometry
		utility::file::file_delete( silent_outfile );
		core::io::silent::BinarySilentStruct pss( opts, ref_pose, "tag" );
		sfd.write_silent_struct( pss, silent_outfile );

		sfd.read_file(silent_outfile);
		TS_ASSERT( sfd.size() > 0 );
		core::io::silent::SilentFileData::iterator iter = sfd.begin();
		iter->fill_pose( restored_pose );
		TS_ASSERT( !restored_pose.is_fullatom() );
		// test rms difference
		Real rms_to_restored = scoring::CA_rmsd( ref_pose, restored_pose );
		TR << "RMS error from centroid save/restore: " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < rms_threshold );

		// test score13 difference
		core::scoring::ScoreFunctionOP scorefxn =
			core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
		Real score_ref = (*scorefxn)(ref_pose);
		Real score_restored = (*scorefxn)(restored_pose);
		Real score_del = std::fabs( score_restored - score_ref );
		TR << "Score difference: " << score_del << std::endl;
		TS_ASSERT( score_del < score_threshold );
		utility::file::file_delete( silent_outfile );
	}

	void test_save_and_restore()
	{
		double rms_threshold = 1e-2;
		double score_threshold = 1e-1;

		pose::Pose ref_pose, restored_pose;
		core::chemical::ResidueTypeSetCOP rsd =
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		core::import_pose::pose_from_file( ref_pose, *rsd, std::string("core/io/bin_silentfile_test.pdb"), core::import_pose::PDB_file);
		std::string const silent_outfile( "core/io/bin_silentfile_test.out" ); // read file w/ non-ideal geometry
		utility::file::file_delete( silent_outfile );
		core::io::silent::SilentFileOptions opts;
		core::io::silent::SilentFileData sfd( opts );
		core::io::silent::BinarySilentStruct pss( opts, ref_pose, "tag" );
		sfd.write_silent_struct( pss, silent_outfile );
		// Read the ProteinSilentStruct from the silent-file

		sfd.read_file( silent_outfile );
		TS_ASSERT( sfd.size() > 0 );
		core::io::silent::SilentFileData::iterator iter = sfd.begin();
		iter->fill_pose( restored_pose, *rsd );

		// test rms difference
		Real rms_to_restored = scoring::CA_rmsd( ref_pose, restored_pose );
		TR << "RMS error from save/restore: " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < rms_threshold );

		// test score13 difference
		core::scoring::ScoreFunctionOP scorefxn =
			core::scoring::ScoreFunctionFactory::create_score_function( "score13_env_hb" );
		Real score_ref = (*scorefxn)(ref_pose);
		Real score_restored = (*scorefxn)(restored_pose);
		Real score_del = std::fabs( score_restored - score_ref );
		TR << "Score difference: " << score_del << std::endl;
		if ( score_del > score_threshold ) {
			restored_pose.dump_pdb( "restored_pose_wrong_score.pdb" );
			ref_pose.dump_pdb( "ref_pose_wrong_score.pdb" );
		}
		TS_ASSERT( score_del < score_threshold );

		// utility::file::file_delete( silent_outfile );

	}

	/// @brief Test whether badly-formatted silent files can be read.
	/// @details The silent file read machinery was recently made more robust to corrupted silent files of the
	/// sort commonly produced by Rosetta for inscrutable reasons.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_read_badly_formatted_silent_file() {
		TR << "Starting test_read_badly_formatted_silent_file." << std::endl;
		TR << "Author: Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)." << std::endl;
		TR << "Failure of this test would indicate that there's a problem reading silent files that have been written by Rosetta with common write issues, like duplicated header information or extra data columns in score lines." << std::endl;

		core::pose::Pose good_pose, questionable_pose;
		core::chemical::ResidueTypeSetCOP rsd( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

		core::io::silent::SilentFileOptions opts;
		core::io::silent::SilentFileData sfd(opts);
		sfd.read_file( "core/io/silent/cyclic_peptide_binary_silent.silent" );
		TS_ASSERT( sfd.size() > 0 );
		TS_ASSERT( sfd.begin()->decoy_tag() == "result_proc0042_34330");
		sfd.begin()->fill_pose( good_pose, *rsd );

		core::io::silent::SilentFileData sfd2(opts);
		sfd2.read_file( "core/io/silent/cyclic_peptide_binary_silent_duplicated_header.silent" );
		TS_ASSERT( sfd2.size() > 0 );
		TS_ASSERT( sfd2.begin()->decoy_tag() == "result_proc0042_34330");
		sfd2.begin()->fill_pose( questionable_pose, *rsd );

		TS_ASSERT_EQUALS( good_pose.size(), 18 );
		TS_ASSERT_EQUALS( good_pose.size(), questionable_pose.size() );

		//std::string const tag1( core::pose::tag_from_pose(good_pose) );
		//std::string const tag2( core::pose::tag_from_pose(questionable_pose) );
		//TR << "TAG1: " << tag1 << "\tTAG2: " << tag2 << "\n";
		//TS_ASSERT( tag1 == tag2 );
		//TS_ASSERT( tag1 == "result_proc0042_34330" );
		//TS_ASSERT( tag2 == "result_proc0042_34330" );

		TR << "\nPHI1\tPHI2\tPSI1\tPSI2\tOMEGA1\tOMEGA2\n";
		for ( core::Size i=1; i<=18; ++i ) {
			TR << good_pose.phi(i) << "\t" << questionable_pose.phi(i) << "\t";
			TR << good_pose.psi(i) << "\t" << questionable_pose.psi(i) << "\t";
			TR << good_pose.omega(i) << "\t" << questionable_pose.omega(i) << "\n";
			TS_ASSERT_DELTA( good_pose.phi( i ), questionable_pose.phi( i ), 0.001 );
			TS_ASSERT_DELTA( good_pose.psi( i ), questionable_pose.psi( i ), 0.001 );
			TS_ASSERT_DELTA( good_pose.omega( i ), questionable_pose.omega( i ), 0.001 );
		}
		TR << "\n";

		core::scoring::ScoreFunctionOP scorefxn( core::scoring::ScoreFunctionFactory::create_score_function( "talaris2014" ) );
		core::Real const score1( (*scorefxn)(good_pose) );
		core::Real const score2( (*scorefxn)(questionable_pose) );
		TR << "SCORE1\tSCORE2\n" << score1 << "\t" << score2 << "\n";
		TS_ASSERT_DELTA( score1, score2, 0.001 );
		TR << std::endl; TR.flush();
	}

	/// @brief Test whether poses that have had their lenghts modified in Rosetta can
	/// be written and read properly.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_save_and_restore_inserting_residues()
	{
		TR << "Starting test_save_and_restore_inserting_residues." << std::endl;
		TR << "Author: Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)." << std::endl;
		TR << "Failure of this test would indicate that there's a problem storing structures whose length Rosetta has modified." << std::endl;
		double rms_threshold = 1e-2;
		double score_threshold = 1e-1;

		pose::Pose ref_pose, restored_pose;
		core::chemical::ResidueTypeSetCOP rsd =
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		core::import_pose::pose_from_file( ref_pose, *rsd, std::string("core/io/bin_silentfile_test.pdb"), core::import_pose::PDB_file);

		core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( rsd->name_map( "ALA" ) ) );
		core::Size nres = ref_pose.size();
		if ( ref_pose.residue(nres).has_variant_type(core::chemical::UPPER_TERMINUS_VARIANT) ) {
			core::pose::remove_variant_type_from_pose_residue( ref_pose, core::chemical::UPPER_TERMINUS_VARIANT, nres );
		}
		ref_pose.append_residue_by_bond(*new_rsd, true, 0, nres, 0, false);
		core::pose::add_variant_type_to_pose_residue( ref_pose, core::chemical::UPPER_TERMINUS_VARIANT, ref_pose.size() );

		std::string const silent_outfile( "core/io/bin_silentfile_test_inserting_residues.out" );
		utility::file::file_delete( silent_outfile );
		core::io::silent::SilentFileOptions opts;
		core::io::silent::SilentFileData sfd(opts);
		core::io::silent::BinarySilentStruct pss( opts, ref_pose, "tag" );
		sfd.write_silent_struct( pss, silent_outfile );
		// Read the ProteinSilentStruct from the silent-file

		sfd.read_file( silent_outfile );
		TS_ASSERT( sfd.size() > 0 );
		core::io::silent::SilentFileData::iterator iter = sfd.begin();
		iter->fill_pose( restored_pose, *rsd );

		// test rms difference
		Real rms_to_restored = scoring::CA_rmsd( ref_pose, restored_pose );
		TR << "RMS error from save/restore: " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < rms_threshold );

		// test score13 difference
		core::scoring::ScoreFunctionOP scorefxn =
			core::scoring::ScoreFunctionFactory::create_score_function( "score13_env_hb" );
		Real score_ref = (*scorefxn)(ref_pose);
		Real score_restored = (*scorefxn)(restored_pose);
		Real score_del = std::fabs( score_restored - score_ref );
		TR << "Score difference: " << score_del << std::endl;
		if ( score_del > score_threshold ) {
			restored_pose.dump_pdb( "restored_pose_wrong_score.pdb" );
			ref_pose.dump_pdb( "ref_pose_wrong_score.pdb" );
		}
		TS_ASSERT( score_del < score_threshold );

		// utility::file::file_delete( silent_outfile );

	}

};
