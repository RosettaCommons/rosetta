// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/io/carbohydrates/pose_io.cxxtest.hh
/// @brief   Test suite for carbohydrate-containing Pose input/output functions.
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/io/carbohydrates/pose_io.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("core.io.carbohydrates.pose_io.cxxtest");

class CarbohydratePoseIOTests : public CxxTest::TestSuite {
public: // Standard methods ///////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init();
	}

	// Destruction
	void tearDown()
	{}


public: // Tests //////////////////////////////////////////////////////////////
	// Confirm that parsing of sugar modifications from carbohydrate code suffixes works properly.
	void test_sugar_modifications_from_suffix()
	{
		using namespace core;
		using namespace io::carbohydrates;

		TR << "Testing sugar_modifications_from_suffix() method." << std::endl;

		utility::vector1< std::pair< core::uint, std::string > > modifications;

		// All of these modifications are taken directly from the list of examples on the IUPAC webpage:
		// www.chem.qmul.ac.uk/iupac/2carb/38.html

		modifications = sugar_modifications_from_suffix( "-ol" );
		TS_ASSERT_EQUALS( modifications.size(), 1 );
		TS_ASSERT_EQUALS( modifications[ 1 ].first, 0 );
		TS_ASSERT_EQUALS( modifications[ 1 ].second, "-ol" );

		modifications = sugar_modifications_from_suffix( "N" );
		TS_ASSERT_EQUALS( modifications.size(), 1 );
		TS_ASSERT_EQUALS( modifications[ 1 ].first, 0 );
		TS_ASSERT_EQUALS( modifications[ 1 ].second, "N" );

		modifications = sugar_modifications_from_suffix( "NAc" );
		TS_ASSERT_EQUALS( modifications.size(), 1 );
		TS_ASSERT_EQUALS( modifications[ 1 ].first, 0 );
		TS_ASSERT_EQUALS( modifications[ 1 ].second, "NAc" );

		modifications = sugar_modifications_from_suffix( "4S" );
		TS_ASSERT_EQUALS( modifications.size(), 1 );
		TS_ASSERT_EQUALS( modifications[ 1 ].first, 4 );
		TS_ASSERT_EQUALS( modifications[ 1 ].second, "S" );

		modifications = sugar_modifications_from_suffix( "N3N" );
		TS_ASSERT_EQUALS( modifications.size(), 2 );
		TS_ASSERT_EQUALS( modifications[ 1 ].first, 0 );
		TS_ASSERT_EQUALS( modifications[ 1 ].second, "N" );
		TS_ASSERT_EQUALS( modifications[ 2 ].first, 3 );
		TS_ASSERT_EQUALS( modifications[ 2 ].second, "N" );

		modifications = sugar_modifications_from_suffix( "A" );
		TS_ASSERT_EQUALS( modifications.size(), 1 );
		TS_ASSERT_EQUALS( modifications[ 1 ].first, 0 );
		TS_ASSERT_EQUALS( modifications[ 1 ].second, "A" );

		modifications = sugar_modifications_from_suffix( "A6Et" );
		TS_ASSERT_EQUALS( modifications.size(), 2 );
		TS_ASSERT_EQUALS( modifications[ 1 ].first, 0 );
		TS_ASSERT_EQUALS( modifications[ 1 ].second, "A" );
		TS_ASSERT_EQUALS( modifications[ 2 ].first, 6 );
		TS_ASSERT_EQUALS( modifications[ 2 ].second, "Et" );

		modifications = sugar_modifications_from_suffix( "5Ac" );
		TS_ASSERT_EQUALS( modifications.size(), 1 );
		TS_ASSERT_EQUALS( modifications[ 1 ].first, 5 );
		TS_ASSERT_EQUALS( modifications[ 1 ].second, "Ac" );

		modifications = sugar_modifications_from_suffix( "2en5Ac" );
		TS_ASSERT_EQUALS( modifications.size(), 2 );
		TS_ASSERT_EQUALS( modifications[ 1 ].first, 2 );
		TS_ASSERT_EQUALS( modifications[ 1 ].second, "en" );
		TS_ASSERT_EQUALS( modifications[ 2 ].first, 5 );
		TS_ASSERT_EQUALS( modifications[ 2 ].second, "Ac" );

		modifications = sugar_modifications_from_suffix( "5Gc" );
		TS_ASSERT_EQUALS( modifications.size(), 1 );
		TS_ASSERT_EQUALS( modifications[ 1 ].first, 5 );
		TS_ASSERT_EQUALS( modifications[ 1 ].second, "Gc" );

		// Technically, this should be 3,4Me2, with a subscripted 2, but I cannot parse subscripts.
		modifications = sugar_modifications_from_suffix( "3,4Me" );
		TS_ASSERT_EQUALS( modifications.size(), 2 );
		TS_ASSERT_EQUALS( modifications[ 1 ].first, 3 );
		TS_ASSERT_EQUALS( modifications[ 1 ].second, "Me" );
		TS_ASSERT_EQUALS( modifications[ 2 ].first, 4 );
		TS_ASSERT_EQUALS( modifications[ 2 ].second, "Me" );

		modifications = sugar_modifications_from_suffix( "2CMe" );
		TS_ASSERT_EQUALS( modifications.size(), 1 );
		TS_ASSERT_EQUALS( modifications[ 1 ].first, 2 );
		TS_ASSERT_EQUALS( modifications[ 1 ].second, "CMe" );

		// Now test for some expected errors.
		TR <<  " -------------- Six (6) error messages about commas should follow. ------------- "  << std::endl;
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( sugar_modifications_from_suffix( "going,ndo" ) );
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( sugar_modifications_from_suffix( "Oxford," ) );
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( sugar_modifications_from_suffix( ",Y" ) );
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( sugar_modifications_from_suffix( "1,Y" ) );
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( sugar_modifications_from_suffix( "X,9" ) );
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( sugar_modifications_from_suffix( "11" ) );
		TR <<  "--------------- Six (6) error messages should have just been output. --------------"  << std::endl;
	}

	// TODO: Add tests for GWS input and output!
};  // class CarbohydratePoseIOTests
