// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/mhc_epitope_energy/MHCEpitopePredictorPreLoaded.cxxtest.hh
/// @brief  Test suite for core::scoring::mhc_epitope_energy::MHCEptiopePredictorPreLoaded, the MHC epitope predictor using a preloaded external SQL or CSV database.
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

// Test headers
#include <cxxtest/TestSuite.h>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorPreLoaded.hh>
#ifndef MULTI_THREADED
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorExternal.hh>
#endif

// Package Headers
#include <test/core/init_util.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("core.scoring.mhc_epitope_energy.MHCEpitopePredictorPreLoaded.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace core;
//using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::mhc_epitope_energy;

class MHCEpitopePredictorPreLoadedTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //
	// Shared variables
	std::string sequence;
	std::string unk_sequence;

	// Shared initialization goes here.
	void setUp() {
		core_init();

		// Setup test sequences
		sequence = "YFCTRAFRILAWIGI";
		unk_sequence = "AAAAAAAAAAAAAAA";
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Check that the MHCEpitopePredictorPreLoaded machinery is working.
	/// @author Brahm Yachnin
	void test_mhc_predictor_preloaded() {
		//Create a new MHCEpitopePredictorPreLoaded object
		MHCEpitopePredictorPreLoadedOP pred_preload( utility::pointer::make_shared<MHCEpitopePredictorPreLoaded>() );

		//Check the default unseen_score_
		TS_ASSERT_EQUALS(pred_preload->get_unseen_score(), 100);

		//Now, load the external database.
		pred_preload->load_database( "core/scoring/mhc_epitope_energy/externaldb.sql" );

		//Check the score for an unknown peptide
		TS_ASSERT_EQUALS(pred_preload->score(unk_sequence), 100);

		//Set the unknown score to something else and re-score
		pred_preload->set_unseen_score(50);
		TS_ASSERT_EQUALS(pred_preload->get_unseen_score(), 50);
		core::Real unk_score = pred_preload->score(unk_sequence);
		TS_ASSERT_EQUALS(unk_score, 50);

		//Check the score of the known sequence
		core::Real known_score = pred_preload->score(sequence);
		TS_ASSERT_EQUALS(known_score, 2);

		//Create a new MHCEpitopePredictorPreLoaded class using a CSV file
		MHCEpitopePredictorPreLoadedOP pred_preload_csv( utility::pointer::make_shared<MHCEpitopePredictorPreLoaded>() );
		pred_preload_csv->load_csv( "core/scoring/mhc_epitope_energy/externaldb.csv" );
		pred_preload_csv->set_unseen_score(50);

		//Check that the CSV scores are the same as the SQL scores
		TS_ASSERT_EQUALS(pred_preload_csv->score(unk_sequence), unk_score);
		TS_ASSERT_EQUALS(pred_preload_csv->score(sequence), known_score);

#ifndef MULTI_THREADED
		//Check if the External predictor gives the same results.
		//Create a new MHCEpitopePredictorExternal object
		MHCEpitopePredictorExternalOP pred_e( utility::pointer::make_shared<MHCEpitopePredictorExternal>() );
		pred_e->connect( "core/scoring/mhc_epitope_energy/externaldb.sql" );
		pred_e->set_unseen_score(50);

		//Check that the CSV scores are the same as the SQL scores
		TS_ASSERT_EQUALS(pred_e->score(unk_sequence), unk_score);
		TS_ASSERT_EQUALS(pred_e->score(sequence), known_score);
#endif

		TR << "End of test_mhc_predictor_preloaded." << std::endl;
	}

};
