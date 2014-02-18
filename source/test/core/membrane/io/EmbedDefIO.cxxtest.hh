// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/io/EmbedDefIO.cxxtest.hh
///
/// @brief 		Test Suite for embedding definitions
/// @details	Tests for embedding definitions!
///
/// @author 	Rebecca Alford

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Tested Classes
#include <core/membrane/io/EmbedDefIO.hh>

#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

#include <protocols/loops/LoopsFileOptions.hh> //this is crazy - lol inspired by Milke :)

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

/// @brief Test Class Embedding Definition Data Loader
class EmbedDefIOTests :  public CxxTest::TestSuite {

public: // tests

	/// @brief Setup Function - Runs at the beginning of each test
	void setUp()
	{
		core_init();

		// IO Class
		edio_ = new core::membrane::io::EmbedDefIO();
	}

	/// @brief Tear Down Function
	void tearDown()
	{}

	/// @brief Test empty file exception
	void testEmptyFile()
	{
		TS_TRACE("Testing empty embed file input");
		TS_ASSERT_THROWS_ANYTHING( edio_->get_embedding_from_file("") );
	}

	/// @brief Testing Expected Parameters
	void testExpected()
	{
		TS_TRACE("Testing expected file initialization");

		// Load data from file
		config_ = edio_->get_embedding_from_file("core/membrane/io/1afo_test.embed");

		// Make Comparisons
		TS_ASSERT_EQUALS( config_->normal_tag.compare("from_topology"), 0 );
		TS_ASSERT_EQUALS( config_->center_tag.compare("from_topology"), 0 );

		TS_ASSERT_EQUALS( config_->normal.x(), -0.0904839 );
		TS_ASSERT_EQUALS( config_->normal.y(), 0.977314 );
		TS_ASSERT_EQUALS( config_->normal.z(), -0.191494 );

		TS_ASSERT_EQUALS( config_->center.x(), 0.98975 );
		TS_ASSERT_EQUALS( config_->center.y(), -0.8255 );
		TS_ASSERT_EQUALS( config_->center.z(), -2.32475 );

		TS_ASSERT_EQUALS( config_->depth, 0 );

	}

private: // data

	// IO Class
	core::membrane::io::EmbedDefIOOP edio_;

	// EmbedConfigInfo Obejct to Load
	core::membrane::util::EmbedConfigInfoOP config_;

};