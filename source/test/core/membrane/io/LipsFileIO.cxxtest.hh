// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file                 core/membrane/io/LipsFileIO.cxxtest.hh
///
/// @brief                 Test Suite for izstream reader class that reads in lips4 file data
/// @details
///
/// @author         Rebecca Alford (rfalford12@gmail.com)

/**
 *  Dependency note - depends upon lips4 file format NOT .lips
 */

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Tested Classes
#include <core/membrane/io/LipoFileIO.hh>

#include <core/membrane/properties/LipidAccInfo.hh>
#include <core/membrane/util/Exceptions.hh>

// Package Headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

/// @brief Test Class: Lipo File loader
class LipsFileIOTests : public CxxTest::TestSuite {

public: // tests

    /// @brief Set Up function - Runs at start of every test
	void setUp()
	{

		using namespace core::membrane::util;
		using namespace core::membrane::io;
        using namespace core::import_pose;

		core_init();
        
		// Create an instance of the lipofile io class
		lfio_ = new core::membrane::io::LipoFileIO();
        
		// Set locator for test lipsfile
		test_lipsfile_ = "core/membrane/io/1afo_test.lips";
	}

	/// @brief Tear Down
	void tearDown()
	{}

	/// @brief Test Empty Input
	void test_emptyLocator()
	{

		TS_TRACE("Testing lips response to empty locator");
		TS_ASSERT_THROWS_ANYTHING( lfio_->get_lips_exp_from_lipofile( "" ) );

	}

	/// @brief Testing invalid file location for lipsfile
	void test_invalidLocator() {
        
		TS_TRACE("Testing lips io throws exception given an invalid locator");
		TS_ASSERT_THROWS_ANYTHING( lfio_->get_lips_exp_from_lipofile( "aaaaa" ) );
	}

	/// @brief Testing Expected test cases - single tmh helix
	void test_correct1afo() {
        
		TS_TRACE("Testing a correctly initialized lips data file");
        
		// Loading a lipid acc file (should use an actual resource description now
		LipidAccInfoOP lips_exp = lfio_->get_lips_exp_from_lipofile( test_lipsfile_ );
        
		// Check for non-null pointer access
		TS_ASSERT( lips_exp );
        
		// Testing quick contents
		TS_ASSERT_EQUALS( lips_exp->lipid_exposure()[26], 0.5);
		TS_ASSERT_EQUALS( lips_exp->lipid_exposure()[15], 0.0);
		TS_ASSERT_EQUALS( lips_exp->lipid_burial()[26], 0.0);
		TS_ASSERT_EQUALS( lips_exp->lipid_burial()[15], 0.0);
	}

private: // data

	// LFIO
	core::membrane::io::LipoFileIOOP lfio_;

	// Test lips file frompose
	std::string test_lipsfile_;
    
};
