// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/mhc_epitope_energy/MHCEpitopePredictorMatrix.cxxtest.hh
/// @brief  Test suite for core::scoring::mhc_epitope_energy::MHCEptiopePredictorMatrix, the MHC epitope predictor using a position weight matrix.
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

// Test headers
#include <cxxtest/TestSuite.h>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorMatrix.hh>

// Package Headers
#include <test/core/init_util.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("core.scoring.mhc_epitope_energy.MHCEpitopePredictorMatrix.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::scoring;
using namespace core::scoring::mhc_epitope_energy;

class MHCEpitopePredictorMatrixTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //
	// Shared variables
	std::string sequence;

	// Shared initialization goes here.
	void setUp() {
		core_init();

		// Setup test sequence
		sequence = "YFCTRAFRI";
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Check that the MHCEpitopePredictorMatrix machinery is working.
	/// @author Brahm Yachnin
	void test_mhc_predictor_matrix() {
		//Create a new MHCEpitopePredictorMatrix object
		MHCEpitopePredictorMatrixOP pred_m( utility::pointer::make_shared<MHCEpitopePredictorMatrix>() );

		//By default, the predictor should have a threshold of 5.0
		TS_ASSERT_EQUALS(pred_m->get_thresh(), 5.0);

		//Now, load the matrix.
		pred_m->load_matrix( std::string("propred8") );

		//Set the threshold to 1.0, and check the score.
		pred_m->set_thresh(1.0);
		TS_ASSERT_EQUALS(pred_m->get_thresh(), 1.0);
		core::Real thresh1_score = pred_m->score(sequence);
		TS_ASSERT_EQUALS( thresh1_score, 1.0 );

		//If that worked, the length of the peptide and the Predictor peptide_length_ should be the same.
		TS_ASSERT_EQUALS( pred_m->get_peptide_length(), sequence.size() );

		//Set the threshold to 10.0, and check the score.
		pred_m->set_thresh(10.0);
		TS_ASSERT_EQUALS(pred_m->get_thresh(), 10.0);
		TS_ASSERT_EQUALS( pred_m->score(sequence), 5.0 );

		//Set the threshold back to 1.0, and check that we get the same score again.
		pred_m->set_thresh(1.0);
		TS_ASSERT_EQUALS(pred_m->get_thresh(), 1.0);
		TS_ASSERT_EQUALS( pred_m->score(sequence), thresh1_score );

		TR << "End of test_mhc_predictor_matrix." << std::endl;
	}

};
