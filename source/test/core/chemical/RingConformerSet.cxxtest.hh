// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	 RingConformerSet.cxxtest.hh
/// @brief   Test suite for ring conformer set building and associated methods
/// @author  Labonte <JWLabonte@jhu.edu>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/chemical/RingConformerSet.hh>

// Utility header
#include <utility/vector1.hh>

// C++ header
#include <string>

class RingConformerSetTests : public CxxTest::TestSuite {
public:
	// Standard methods ///////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		using namespace std;
		using namespace utility;
		using namespace core::chemical;

		core_init_with_additional_options( "-out:levels core.chemical.RingConformerSet:400" );

		vector1< string > lows;

		lows.push_back( "1E" );
		lows.push_back( "2E" );
		lows.push_back( "3E" );
		lows.push_back( "PINEAPPLE" );  // should be ignored by RingConformerSet

		set5_ = core::chemical::RingConformerSetOP( new RingConformerSet( 5, "1E", lows ) );
		set6_ = core::chemical::RingConformerSetOP( new RingConformerSet( 6, "", lows ) );  // should default to 4C1
	}

	// Destruction
	void tearDown()
	{}

	// Tests //////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Confirm that the proper number of conformers are loaded into 5- and 6-membered ring sets.
	void test_get_all_nondegenerate_conformers()
	{
		TS_TRACE(
				"Testing get_all_nondegenerate_conformers() method of RingConformerSet for 5- and 6-membered rings." );
		TS_ASSERT_EQUALS( set5_->get_all_nondegenerate_conformers().size(), 20 );
		TS_ASSERT_EQUALS( set6_->size(), 38 );
	}

	// Confirm that the proper conformers are loaded by get_ideal_conformer_by_name() and
	// get_ideal_conformer_by_CP_parameters.
	void test_get_ideal_conformer_by_descriptor_methods()
	{
		using namespace core;
		using namespace core::chemical;
		using namespace utility;

		TS_TRACE( "Testing get_ideal_conformer_by_name() and get_ideal_conformer_by_CP_parameters() methods of"
				" RingConformerSet for 5- and 6-membered rings." );

		vector1< Real > params5, params6, params_bad;
		params5.resize( 2 );
		params6.resize( 3 );

		params5[ q ] = 0.4;
		params5[ PHI ] = 180.0;

		params6[ q ] = 0.55;
		params6[ PHI ] = 180.0;
		params6[ THETA ] = 90.0;

		TS_ASSERT_EQUALS( set5_->get_ideal_conformer_by_name( "EO" ).specific_name,
				set5_->get_ideal_conformer_by_CP_parameters( params5 ).specific_name );
		TS_ASSERT_EQUALS( set6_->get_ideal_conformer_by_name( "BO,3" ).specific_name,
				set6_->get_ideal_conformer_by_CP_parameters( params6 ).specific_name );

		// Test that rounding is handled properly.
		params5[ PHI ] = 43.2;  // should round to 36.0

		TS_ASSERT_EQUALS( set5_->get_ideal_conformer_by_name( "E1" ).specific_name,
				set5_->get_ideal_conformer_by_CP_parameters( params5 ).specific_name );

		params6[ PHI ] = 123.4;  // should round to 120.0
		params6[ THETA ] = 123.4;  // should round to 135.0

		TS_ASSERT_EQUALS( set6_->get_ideal_conformer_by_name( "5E" ).specific_name,
				set6_->get_ideal_conformer_by_CP_parameters( params6 ).specific_name );

		// Test that chairs are handled correctly, as they reside at the poles where phi is meaningless.
		params6[ THETA ] = 0.0;

		TS_ASSERT_EQUALS( set6_->get_ideal_conformer_by_name( "4C1" ).specific_name,
				set6_->get_ideal_conformer_by_CP_parameters( params6 ).specific_name );

		// Test for bad input.
		TS_TRACE( "A not-found error should follow:" );
		try {
			set5_->get_ideal_conformer_by_name( "FOO" );  // This should force an exit.
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_Base const & e) {
			TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ) ),
					"ERROR: No conformer with given name found in this set; exiting.\n\n" );
			TS_TRACE( "The above error message was expected." );
		}
		TS_TRACE( "A series of bad-input errors should follow:" );
		try {
			set6_->get_ideal_conformer_by_CP_parameters( params_bad );  // This should force an exit.
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_Base const & e) {
			TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ) ),
					"ERROR: An N-membered ring is described by exactly N-3 Cremer-Pople parameters, "
					"yet a different number was provided; exiting.\n\n" );
			TS_TRACE( "The above error message was expected." );
		}
		params_bad.push_back( 0.0 );  // still bad because not enough params for a 6-membered ring
		try {
			set6_->get_ideal_conformer_by_CP_parameters( params_bad );  // This should force an exit.
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_Base const & e) {
			TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ) ),
					"ERROR: An N-membered ring is described by exactly N-3 Cremer-Pople parameters, "
					"yet a different number was provided; exiting.\n\n" );
			TS_TRACE( "The above error message was expected." );
		}
		params_bad.push_back( 180.0 );
		params_bad.push_back( 90.0 );  // still bad because q = 0.0 (planar)
		try {
			set6_->get_ideal_conformer_by_CP_parameters( params_bad );  // This should force an exit.
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_Base const & e) {
			TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ) ),
					"ERROR: Planar ring conformations are not handled by Rosetta; "
					"please specify a non-zero q value; exiting.\n\n" );
			TS_TRACE( "The above error message was expected." );
		}
		params_bad.push_back( 90.0 );  // still bad because too many params
		try {
			set6_->get_ideal_conformer_by_CP_parameters( params_bad );  // This should force an exit.
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_Base const & e) {
			TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ) ),
					"ERROR: An N-membered ring is described by exactly N-3 Cremer-Pople parameters, "
					"yet a different number was provided; exiting.\n\n" );
			TS_TRACE( "The above error message was expected." );
		}
	}

	// Confirm that the sub-sets have been loaded properly.
	// Note: This tests a method that returns a random value from a subset; as such, it is possible for it to pass on
	// occasion when it should fail.  For example, if "PINEAPPLE" were incorrectly assigned to the subset, this could
	// go undetected, but not for long, so I don't think this is a real problem.
	void test_get_conformer_methods()
	{
		using namespace std;
		using namespace utility;

		TS_TRACE( "Testing that RingConformer subsets have been populated correctly." );
		TS_ASSERT_EQUALS( set5_->get_lowest_energy_conformer().specific_name, "1E" );
		TS_ASSERT_EQUALS( set6_->get_lowest_energy_conformer().specific_name, "4C1" );

		vector1< string > expected_lows;
		expected_lows.push_back( "1E" );
		expected_lows.push_back( "2E" );
		expected_lows.push_back( "3E" );

		core::chemical::RingConformer random_conformer_from_low_energy_subset;

		random_conformer_from_low_energy_subset = set5_->get_random_local_min_conformer();
		TS_ASSERT( expected_lows.contains( random_conformer_from_low_energy_subset.specific_name ) );

		expected_lows.push_back( "4C1" );
		random_conformer_from_low_energy_subset = set6_->get_random_local_min_conformer();
		TS_ASSERT( expected_lows.contains( random_conformer_from_low_energy_subset.specific_name ) );
	}

private:
	// Private data ////////////////////////////////////////////////////////////
	core::chemical::RingConformerSetOP set5_;
	core::chemical::RingConformerSetOP set6_;

};  // class RingConformerSetTests
