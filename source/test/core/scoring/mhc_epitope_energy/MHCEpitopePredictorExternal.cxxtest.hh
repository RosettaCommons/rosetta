// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/mhc_epitope_energy/MHCEpitopePredictorExternal.cxxtest.hh
/// @brief  Test suite for core::scoring::mhc_epitope_energy::MHCEptiopePredictorExternal, the MHC epitope predictor using an external SQL database.
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

// Test headers
#include <cxxtest/TestSuite.h>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorExternal.hh>

// Package Headers
#include <test/core/init_util.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("core.scoring.mhc_epitope_energy.MHCEpitopePredictorExternal.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace core;
//using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::mhc_epitope_energy;

class MHCEpitopePredictorExternalTests : public CxxTest::TestSuite {

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

	/// @brief Check that the MHCEpitopePredictorExternal machinery is working.
	/// @author Brahm Yachnin
	void test_mhc_predictor_external() {
#ifndef MULTI_THREADED
		//Create a new MHCEpitopePredictorExternal object
		MHCEpitopePredictorExternalOP pred_e( utility::pointer::make_shared<MHCEpitopePredictorExternal>() );

		//Check the default unseen_penalty_
		TS_ASSERT_EQUALS(pred_e->get_unseen_score(), 100);

		//Check that there is no database
		TS_ASSERT(! pred_e->get_database() );

		//Now, load the external database.
		pred_e->connect( "core/scoring/mhc_epitope_energy/externaldb.sql" );
		TS_ASSERT( pred_e->get_database() );

		//Check the score for an unknown peptide
		TS_ASSERT_EQUALS(pred_e->score(unk_sequence), 100);

		//Set the unknown score to something else and re-score
		pred_e->set_unseen_score(50);
		TS_ASSERT_EQUALS(pred_e->get_unseen_score(), 50);
		TS_ASSERT_EQUALS(pred_e->score(unk_sequence), 50);

		//Check the score of the known sequence
		TS_ASSERT_EQUALS(pred_e->score(sequence), 2);

#endif
		TR << "End of test_mhc_predictor_external." << std::endl;
	}

};
