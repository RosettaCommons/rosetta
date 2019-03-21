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

// Package Headers
#include <test/core/init_util.hh>

#include <basic/Tracer.hh>

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
};
