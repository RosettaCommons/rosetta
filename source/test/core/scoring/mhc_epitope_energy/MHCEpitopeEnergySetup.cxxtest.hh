// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/mhc_epitope_energy/MHCEpitopeEnergy.cxxtest.hh
/// @brief  Test suite for core::scoring::mhc_epitope_energy::MHCEptiopeEnergySetup, encapsulating prediction of MHC-peptide binding, to enable deimmunization by mutagenic epitope deletion
/// The code is largely based on (via copying and modifying) NetChargeEnergy (helpers) and HBNetEnergy (only updating around substitution position).
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

// Test headers
#include <cxxtest/TestSuite.h>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergySetup.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorPreLoaded.hh>

// Package Headers
#include <test/core/init_util.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <utility/io/izstream.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>


static basic::Tracer TR("core.scoring.mhc_epitope_energy.MHCEpitopeEnergySetup.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::scoring;

using namespace core::scoring::mhc_epitope_energy;

class MHCEpitopeEnergySetupTests_initialization : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Parse the mhc_epitope parameters from an .mhc file format, and test them.
	/// @brief This uses a matrix-type setup
	/// @author Brahm Yachnin
	void test_mhc_init_matrix( ) {
		// Make a new MHCEpitopeEnergySetup object.
		MHCEpitopeEnergySetupOP mhc_from_file( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );

		// The setup object should be the default, always-zero predictor.
		TS_ASSERT( mhc_from_file->is_default() );
		TS_ASSERT_EQUALS( mhc_from_file->raw_score("FRILAAIGI"), 0);

		// Initialize from a .mhc file
		mhc_from_file->initialize_from_file( "propred8_5.mhc" );

		// The setup object should no longer be the default.
		TS_ASSERT( ! mhc_from_file->is_default() );

		// We are using propred, so the peptide length should be 9.
		TS_ASSERT_EQUALS( mhc_from_file->get_peptide_length(), 9 );
		// Get the raw score from an arbitrary peptide using the predictor.
		TS_ASSERT_EQUALS( mhc_from_file->raw_score("FRILAAIGI"), 8.0 );

		// Make another MHCEpitopeEnergySetup object, this time to config using "file contents."
		MHCEpitopeEnergySetupOP mhc_from_string( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );

		// Set up a string object that matches the .mhc file
		std::string config_string = "method matrix propred8\nalleles *\nthreshold 5";

		// Initialize it using the string
		mhc_from_string->initialize_from_file_contents( config_string );

		// The file and file contents predictors should be the same.
		TS_ASSERT_EQUALS( mhc_from_file->report(), mhc_from_string->report() ); //I think this is the best we can do without comparing all accessors.

		// Test using operator==
		TS_ASSERT_EQUALS ( *mhc_from_file == *mhc_from_string, true );

		// One with a different threshold
		MHCEpitopeEnergySetupOP mhc_from_string2( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );
		std::string config_string2 = "method matrix propred8\nalleles *\nthreshold 2";
		mhc_from_string2->initialize_from_file_contents( config_string2 );
		TS_ASSERT_EQUALS ( *mhc_from_string2 == *mhc_from_string, false );

		// TODO: when have ability to change set of alleles, test that contribution to operator==

		// No predictor
		MHCEpitopeEnergySetupOP mhc_from_string3( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );
		TS_ASSERT_EQUALS ( *mhc_from_string3 == *mhc_from_string, false );

		TR << "End of test_mhc_init_matrix." << std::endl;
	}

	/// @brief Parse the mhc_epitope parameters from an .mhc file format, and test them.
	/// @brief This uses a external database-type setup
	/// @author Brahm Yachnin
	void test_mhc_init_database( ) {
#ifndef MULTI_THREADED
		// Make a new MHCEpitopeEnergySetup object.
		MHCEpitopeEnergySetupOP mhc_from_file( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );

		// The setup object should be the default, always-zero predictor.
		TS_ASSERT( mhc_from_file->is_default() );

		// Initialize from a .mhc file
		mhc_from_file->initialize_from_file( "core/scoring/mhc_epitope_energy/external_db.mhc" );

		// The setup object should no longer be the default.
		TS_ASSERT( ! mhc_from_file->is_default() );

		// We are using a NetMHCII database, so the peptide length should be 15.
		TS_ASSERT_EQUALS( mhc_from_file->get_peptide_length(), 15 );
		// Get the raw score from an arbitrary peptide using the predictor.
		TS_ASSERT_EQUALS( mhc_from_file->raw_score("YFCTRAFRILAWIGI"), 2.0 );

		// Make another MHCEpitopeEnergySetup object, this time to config using "file contents."
		MHCEpitopeEnergySetupOP mhc_from_string( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );

		// Set up a string object that matches the .mhc file
		std::string config_string = "method external core/scoring/mhc_epitope_energy/externaldb.sql\nalleles *\nunseen penalize 100";

		// Initialize it using the string
		mhc_from_string->initialize_from_file_contents( config_string );

		// The file and file contents predictors should be the same.
		TS_ASSERT_EQUALS( mhc_from_file->report(), mhc_from_string->report() ); //I think this is the best we can do without comparing all accessors.

		// Test using operator==
		TS_ASSERT_EQUALS ( *mhc_from_file == *mhc_from_string, true );

		// One with a different unseen penalty
		MHCEpitopeEnergySetupOP mhc_from_string2( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );
		std::string config_string2 = "method external core/scoring/mhc_epitope_energy/externaldb.sql\nalleles *\nunseen penalize 1";
		mhc_from_string2->initialize_from_file_contents( config_string2 );
		TS_ASSERT_EQUALS ( *mhc_from_string2 == *mhc_from_string, false );

		// Different type of predictor
		MHCEpitopeEnergySetupOP mhc_from_string3( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );
		std::string config_string3 = "method matrix propred8\nalleles *\nthreshold 2";
		mhc_from_string3->initialize_from_file_contents( config_string3 );
		TS_ASSERT_EQUALS ( *mhc_from_string3 == *mhc_from_string, false );

#endif
		TR << "End of test_mhc_init_database." << std::endl;
	}

	/// @brief Parse the mhc_epitope parameters from an .mhc file format, and test them.
	/// @brief This uses a preloaded csv-type setup
	/// @author Brahm Yachnin
	void test_mhc_init_csv( ) {
		// Make a new MHCEpitopeEnergySetup object.
		MHCEpitopeEnergySetupOP mhc_from_file( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );

		// The setup object should be the default, always-zero predictor.
		TS_ASSERT( mhc_from_file->is_default() );

		// Initialize from a .mhc file
		mhc_from_file->initialize_from_file( "core/scoring/mhc_epitope_energy/preload_csv.mhc" );

		// The setup object should no longer be the default.
		TS_ASSERT( ! mhc_from_file->is_default() );

		// We are using a NetMHCII database, so the peptide length should be 15.
		TS_ASSERT_EQUALS( mhc_from_file->get_peptide_length(), 15 );
		// Get the raw score from an arbitrary peptide using the predictor.
		TS_ASSERT_EQUALS( mhc_from_file->raw_score("YFCTRAFRILAWIGI"), 2.0 );

		// Make another MHCEpitopeEnergySetup object, this time to config using "file contents."
		MHCEpitopeEnergySetupOP mhc_from_string( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );

		// Set up a string object that matches the .mhc file
		std::string config_string = "method preloaded csv core/scoring/mhc_epitope_energy/externaldb.csv\nalleles *\nunseen score 100";

		// Initialize it using the string
		mhc_from_string->initialize_from_file_contents( config_string );

		// The file and file contents predictors should be the same.
		TS_ASSERT_EQUALS( mhc_from_file->report(), mhc_from_string->report() ); //I think this is the best we can do without comparing all accessors.

		// Test using operator==
		TS_ASSERT_EQUALS ( *mhc_from_file == *mhc_from_string, true );

		// One with a different unseen score
		MHCEpitopeEnergySetupOP mhc_from_string2( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );
		std::string config_string2 = "method preloaded csv core/scoring/mhc_epitope_energy/externaldb.csv\nalleles *\nunseen penalize 1";
		mhc_from_string2->initialize_from_file_contents( config_string2 );
		TS_ASSERT_EQUALS ( *mhc_from_string2 == *mhc_from_string, false );

		// Set up a predictor that ignores unseen peptides
		MHCEpitopeEnergySetupOP mhc_ignore_unseen( utility::pointer::make_shared< MHCEpitopeEnergySetup >() );
		std::string config_ignore_unseen = "method preloaded csv core/scoring/mhc_epitope_energy/externaldb.csv\nalleles *\nunseen ignore";
		mhc_ignore_unseen->initialize_from_file_contents( config_ignore_unseen );
		TS_ASSERT_EQUALS( mhc_ignore_unseen->raw_score( "AFRILAWIGIQNPTS" ), 4 ); // A peptide that is known should have a non-zero score.
		TS_ASSERT_EQUALS( mhc_ignore_unseen->raw_score( "AAAAAAWIGIQNPTS" ), 0 ); // A peptide that is UNknown should have a zero score.

		TR << "End of test_mhc_init_csv." << std::endl;
	}

	/// @brief Use the raw xform settings and test that raw_score and xform work.
	/// @author Brahm Yachnin
	void test_mhc_xformed_raw( ) {
		// Initialize another MHCEpitopeEnergySetup object, with a raw xform.
		core::Real raw_offset = 5.0;
		std::string xform_raw = "method matrix propred8\nalleles *\nthreshold 5\nxform raw " + std::to_string(raw_offset);
		MHCEpitopeEnergySetupOP mhc_raw_xform( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );
		mhc_raw_xform->initialize_from_file_contents( xform_raw );

		// Calculate the raw score with this predictor, which doesn't include the offset.
		core::Real rawscore = mhc_raw_xform->raw_score("FRILAAIGI");
		// Calculate the xformed score, and make sure it is equal to rawscore - raw_offset.  (Using rawscore as native too, as it doesn't matter here.)
		core::Real xformed_raw = mhc_raw_xform->xform( rawscore, rawscore );
		TS_ASSERT_EQUALS( xformed_raw, rawscore - raw_offset );
		// Do the same thing using score in one step.
		TS_ASSERT_EQUALS( xformed_raw, mhc_raw_xform->score( "FRILAAIGI", rawscore ) );

		// Repeat the raw xform, but with a larger offset that exceeds the score.  xform should return 0.
		raw_offset = 10.0;
		xform_raw = "method matrix propred8\nalleles *\nthreshold 5\nxform raw " + std::to_string(raw_offset);
		mhc_raw_xform = utility::pointer::make_shared< MHCEpitopeEnergySetup >();
		mhc_raw_xform->initialize_from_file_contents( xform_raw );

		// Calculate the raw score with this predictor, which doesn't include the offset.
		rawscore = mhc_raw_xform->raw_score("FRILAAIGI");
		// Calculate the xformed score, and make sure it returns 0.  (Using rawscore as native too, as it doesn't matter here.)
		core::Real xformed_raw_zeroed = mhc_raw_xform->xform( rawscore, rawscore );
		TS_ASSERT_EQUALS( xformed_raw_zeroed, 0 );
		// Do the same thing using score in one step.
		TS_ASSERT_EQUALS( xformed_raw_zeroed, mhc_raw_xform->score( "FRILAAIGI", rawscore ) );

		TR << "End of test_mhc_xformed_raw." << std::endl;
	}

	/// @brief Use the relative+ xform settings and test that raw_score and xform work.
	/// @author Brahm Yachnin
	void test_mhc_xformed_relative_additive( ) {
		// Set a value for native to be 6.
		core::Real native = 6;

		// Initialize another MHCEpitopeEnergySetup object, with a relative+ xform.
		core::Real rel_offset = 1.0;
		std::string xform_rel = "method matrix propred8\nalleles *\nthreshold 5\nxform relative+ " + std::to_string(rel_offset);
		MHCEpitopeEnergySetupOP mhc_rel_xform( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );
		mhc_rel_xform->initialize_from_file_contents( xform_rel );

		// Calculate the raw score with this predictor, which doesn't include the offset.
		core::Real rawscore = mhc_rel_xform->raw_score("FRILAAIGI");
		// Calculate the xformed score, and make sure it is equal to rawscore - native + rel_offset.
		core::Real xformed_rel = mhc_rel_xform->xform( rawscore, native );
		TS_ASSERT_EQUALS( xformed_rel, rawscore - native + rel_offset );
		// Do the same thing using score in one step.
		TS_ASSERT_EQUALS( xformed_rel, mhc_rel_xform->score( "FRILAAIGI", native ) );

		// Repeat the rel xform, but changing the native such that it is greater than current.  xform should return 0.
		native = 10;

		// Calculate the xformed score, and make sure it returns 0.
		core::Real xformed_rel_zeroed = mhc_rel_xform->xform( rawscore, native );
		TS_ASSERT_EQUALS( xformed_rel_zeroed, 0 );
		// Do the same thing using score in one step.
		TS_ASSERT_EQUALS( xformed_rel_zeroed, mhc_rel_xform->score( "FRILAAIGI", native ) );

		TR << "End of test_mhc_xformed_relative_additive." << std::endl;
	}

	/// @brief Use the relative* xform settings and test that raw_score and xform work.
	/// @author Brahm Yachnin
	void test_mhc_xformed_relative_multiplicative( ) {
		// Set a value for native to be 6.
		core::Real native = 6;

		// Initialize another MHCEpitopeEnergySetup object, with a relative* xform.
		core::Real rel_offset = 1.1;
		std::string xform_rel = "method matrix propred8\nalleles *\nthreshold 5\nxform relative* " + std::to_string(rel_offset);
		MHCEpitopeEnergySetupOP mhc_rel_xform( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );
		mhc_rel_xform->initialize_from_file_contents( xform_rel );

		// Calculate the raw score with this predictor, which doesn't include the offset.
		core::Real rawscore = mhc_rel_xform->raw_score("FRILAAIGI");
		// Calculate the xformed score, and make sure it is equal to rawscore - native * rel_offset.
		core::Real xformed_rel = mhc_rel_xform->xform( rawscore, native );
		TS_ASSERT_EQUALS( xformed_rel, rawscore - native * rel_offset );
		// Do the same thing using score in one step.
		TS_ASSERT_EQUALS( xformed_rel, mhc_rel_xform->score( "FRILAAIGI", native ) );

		// Repeat the rel xform, but changing the native such that it is greater than current.  xform should return 0.
		native = 7.9;

		// Calculate the xformed score, and make sure it returns 0.
		core::Real xformed_rel_zeroed = mhc_rel_xform->xform( rawscore, native );
		TS_ASSERT_EQUALS( xformed_rel_zeroed, 0 );
		// Do the same thing using score in one step.
		TS_ASSERT_EQUALS( xformed_rel_zeroed, mhc_rel_xform->score( "FRILAAIGI", native ) );

		TR << "End of test_mhc_xformed_relative_multiplicative." << std::endl;
	}

	/// @brief This test checks that the iedb_data.db file in the database is valid.
	/// @author Brahm Yachnin
	/// @details This makes some "spot-checks" on database entries, but updates to IEDB may make some of these fail.
	/// @details Updating the peptide list to account for underlying changes in the database is acceptable, as long as database integrity is carefully checked.
	void test_iedb_database_integrity() {
		//Create a new MHCEpitopeEnergySetup class object
		MHCEpitopeEnergySetupOP iedb_predictor( utility::pointer::make_shared<MHCEpitopeEnergySetup>() );

		//Initialize from the database iedb_user.mhc config file
		iedb_predictor->initialize_from_file("iedb_data.mhc");

		//Grab the scoring map from the Predictor
		MHCEpitopePredictorPreLoadedOP iedb_predictor_pred = utility::pointer::dynamic_pointer_cast<MHCEpitopePredictorPreLoaded>(iedb_predictor->get_predictor());
		TS_ASSERT(iedb_predictor_pred);
		std::map< std::string, core::Real > scores_map = iedb_predictor_pred->get_scores_map();

		//The IEDB database probably should only get bigger, not smaller.  As such, this map should only get bigger with database updates.
		//Add a unit test to make sure that the database is at least as big as it was when downloaded on June 21, 2019.
		//As the database grows, this test can be incremented accordingly.  A smaller size should be double-checked carefully.
		const core::Size min_db_size = 190594; //Can be increased, but probably not decreased.
		TR << "Current iedb_data.db database size: " << scores_map.size() << std::endl;
		TS_ASSERT_LESS_THAN_EQUALS(min_db_size, scores_map.size() );

		//Let's check a few peptides that are in the database.
		//Because we have set up the IEDB database to score binders as 1 and non-binders as 0, all binders in the db should be 1.

		//First, set up vectors of assay_mhc_ligand_binding peptides and assay_mhc_ligand_elution peptides (15mers).
		//Additional peptides in the database can be added here, though there isn't a compelling reason to do so.
		const utility::vector1 < std::string > binding_peps = {"FYKETSYTT", "SHCNEMSWI", "YLVSTQEFR", "LLVLAVLCS", "FVSPTPGQR", "APGGAKKPL", "SWQTYVDEH", "KAFLETSNN", "TSTRTYSLG", "MKRYSAPSE", "FMLVAHYAI", "EDLGNCLNK", "QPSSGNYVH", "AGCSEQEVN", "LYAVATTIL"};
		const utility::vector1 < std::string > elution_peps = {"NSSPSAKDI", "AKLSYYDGM", "GPGVPQASG", "EEHASADVE", "DMKNVPEAF", "LNSIKDVEQ", "DLQTIHSRE", "YSLSSVVTV", "YRLHRALDA", "KEDGVITAS"};
		//Finally, set up the same for unseen peptides.  In theory, some of these peptides could appear in the db and make the test fail.
		//As long as the peptide has really been added to the database, this is OK.
		const utility::vector1 < std::string > unknown_peps = {"FYAETSYWT", "SHWNEDSWI", "YEVSTNEFR", "LDVLAVNCS", "FVWPTPMQR", "DCKNVPSAF", "LQSIKDVDQ", "DLYTIHTRE", "YSISSVWTV", "YRLRRALDM"};

		//Starting off with assay_mhc_ligand_binding peptides
		for ( core::Size i=1; i <= binding_peps.size(); ++i ) {
			TS_ASSERT_EQUALS(iedb_predictor->raw_score(binding_peps[i]), 1);
		}

		//Now, the assay_mhc_ligand_elution peptides
		for ( core::Size i=1; i <= elution_peps.size(); ++i ) {
			TS_ASSERT_EQUALS(iedb_predictor->raw_score(elution_peps[i]), 1);
		}

		//Now, some peptides that shouldn't be in the database.  These should have scores of 0.
		for ( core::Size i=1; i <= unknown_peps.size(); ++i ) {
			TS_ASSERT_EQUALS(iedb_predictor->raw_score(unknown_peps[i]), 0);
		}

		//Test that the .info file is OK.  We will check if the file exists, has the right type of info.
		//It is still possible for a user to commit only the db file, and not the info file.
		//This can make this test seem to pass, even though it shouldn't.
		const std::string info_file = basic::database::full_name("scoring/score_functions/mhc_epitope/iedb_data.info");
		const std::string db_file = basic::database::full_name("scoring/score_functions/mhc_epitope/iedb_data.db");

		//Check that iedb_data.info exists.
		utility::io::izstream info_filestream;
		info_filestream.open(info_file);
		TS_ASSERT(info_filestream);

		//Check that the file has text consistent with a successful database refresh.
		//We first check that the first line is not empty, an that it does NOT match the first line of an 'in-progress' db refresh.
		//Since the first line of a successfully refreshed db .info file varies (git username and date), check for a match with the second line.
		std::string firstline = "";
		std::string secondline = "";
		utility::io::getline(info_filestream, firstline); //Get the first line of the .info file.
		utility::io::getline(info_filestream, secondline); //Get the second line of the .info file.
		TS_ASSERT_DIFFERS(firstline, "");
		TS_ASSERT_DIFFERS(firstline, "DO NOT COMMIT THIS FILE AND iedb_data.db.  There may be a problem with iedb_data.db.");
		TS_ASSERT_EQUALS(secondline, "The database was then used to generate iedb_data.db using the following command:");

		info_filestream.close();

		TR << "End of test_iedb_database_integrity." << std::endl;
	}
};
