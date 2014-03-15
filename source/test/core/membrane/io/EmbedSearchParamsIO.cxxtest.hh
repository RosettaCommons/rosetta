// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/io/EmbedSearchParamsIO.cxxtest.hh
///
/// @brief 		Test Suite for loading an embedding search params object form options
/// @details
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

/**
 * Pretty simple test - should just make sure the object initializes correctly
 * from options
 */

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Tested Classes
#include <core/membrane/io/EmbedSearchParamsOptions.hh>
#include <core/membrane/io/EmbedSearchParamsIO.hh>
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/definitions_util.hh>

// Package Headers
#include <core/types.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

class EmbedSearchParamsIOTests : public CxxTest::TestSuite {
public: // tests

	/// @brief Setup - runs before each test
	void setUp()
	{
		// Initialize Test Suites
		core_init();
	}

	/// @brief Tear Down - runs after each test
	void tearDown()
	{}

	/// @brief test correct initialization
	void test_correctInit()
	{
		TS_TRACE("Testing correct initialization of embedding search parameters");

        esp_io_ = new core::membrane::io::EmbedSearchParamsIO();
		esp_opts_ = new core::membrane::io::EmbedSearchParamsOptions();
        
		// Initializing esp opts
		esp_opts_->set_normal_search(true);
		esp_opts_->set_normal_start_angle(10);
		esp_opts_->set_normal_delta_angle(10);
		esp_opts_->set_normal_max_angle(10);
		esp_opts_->set_center_search(true);
		esp_opts_->set_center_max_delta(10);
		esp_opts_->set_normal_mag(20);
		esp_opts_->set_center_mag(30);
		esp_opts_->set_normal_cycles(10);
        esp_opts_->set_penalties(true);
        esp_opts_->set_no_interpolate_mpair(true);
        
        // Create opts ref obj
        core::membrane::io::EmbedSearchParamsOptions const & opts ( *esp_opts_ );
        
        // Get params from file
        core::membrane::util::EmbedSearchParamsOP params = esp_io_->get_embed_params_from_file(opts);

		// Checking Initialization
		TS_ASSERT_EQUALS( params->normal_search, true );
		TS_ASSERT_EQUALS( params->normal_start_angle, 10 );
		TS_ASSERT_EQUALS( params->normal_delta_angle, 10 );
		TS_ASSERT_EQUALS( params->normal_max_angle, 10 );
		TS_ASSERT_EQUALS( params->center_search, true );
		TS_ASSERT_EQUALS( params->center_max_delta, 10 );
		TS_ASSERT_EQUALS( params->normal_mag, 20 );
		TS_ASSERT_EQUALS( params->center_mag, 30 );
		TS_ASSERT_EQUALS( params->normal_cycles, 10 );
        TS_ASSERT_EQUALS( params->penalties, true);
        TS_ASSERT_EQUALS( params->no_interpolate_Mpair, true);
	}

private: // test data

	// IO Class
	core::membrane::io::EmbedSearchParamsIOOP esp_io_;

	// Correct Opts Pointer
	core::membrane::io::EmbedSearchParamsOptionsOP esp_opts_;

};
