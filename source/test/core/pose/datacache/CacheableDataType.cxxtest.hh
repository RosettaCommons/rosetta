// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/datacache/CachableDataType.cxxtest.hh
/// @brief  test suite for core::pose::datacache::CachableDataType
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package headers
#include <core/pose/datacache/CacheableDataType.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("core.pose.datacache.CachableDataType.cxxtest");

// --------------- Test Class --------------- //

class CachableDataTypeTests : public CxxTest::TestSuite {

public:

	// shared initialization
	void setUp() {
	}

	// shared finalization
	void tearDown() {
	}

	// --------------- Test Cases --------------- //

	/// @brief Test the get name function.
	/// This also serves to exercise the runtime checks within the name map intializer
	void test_get_name() {
		using namespace core::pose::datacache;

		// Spot check some values.
		TS_ASSERT_EQUALS( CacheableDataType::get_name( CacheableDataType::ARBITRARY_FLOAT_DATA ), "ARBITRARY_FLOAT_DATA");
		TS_ASSERT_EQUALS( CacheableDataType::get_name( CacheableDataType::MEMBRANE_EMBED ), "MEMBRANE_EMBED");
		TS_ASSERT_EQUALS( CacheableDataType::get_name( CacheableDataType::FULL_MODEL_INFO ), "FULL_MODEL_INFO");
		TS_ASSERT_EQUALS( CacheableDataType::get_name( CacheableDataType::NMR_RDC_DATA ), "NMR_RDC_DATA");

		// Now check that we can get all of them without erroring out and without getting an empty string.
		for ( int ii(1); ii <= CacheableDataType::num_cacheable_data_types; ++ii ) {
			TS_ASSERT_DIFFERS( CacheableDataType::get_name( CacheableDataType::Enum(ii) ), "" );
		}
	}

};
