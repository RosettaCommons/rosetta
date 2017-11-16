// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/chemical/rings/RingConformerSet.cxxtest.hh
/// @brief   Test suite for ring conformer set building and associated methods
/// @author  Labonte <JWLabonte@jhu.edu>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/chemical/rings/RingConformerSet.hh>

// Utility header
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ header
#include <string>

static basic::Tracer TR("core.chemical.rings.RingConformerSet.cxxtest");

class RingConformerSetTests : public CxxTest::TestSuite {
public:
	// Standard methods ///////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		using namespace std;
		using namespace utility;
		using namespace core::chemical;

		core_init_with_additional_options( "-out:levels core.chemical.rings.RingConformerSet:400" );

		vector1< string > lows;

		lows.push_back( "1E" );
		lows.push_back( "2E" );
		lows.push_back( "3E" );
		lows.push_back( "PINEAPPLE" );  // should be ignored by RingConformerSet

		set5_ = rings::RingConformerSetOP( new rings::RingConformerSet( 5, "1E", lows ) );
		set6_ = rings::RingConformerSetOP( new rings::RingConformerSet( 6, "", lows ) );  // should default to 4C1
	}

	// Destruction
	void tearDown()
	{}

	// Tests //////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Confirm that the proper number of conformers are loaded into 5- and 6-membered ring sets.
	void test_get_all_nondegenerate_conformers()
	{
		TR <<
			"Testing get_all_nondegenerate_conformers() method of RingConformerSet for 5- and 6-membered rings." << std::endl;
		TS_ASSERT_EQUALS( set5_->get_all_nondegenerate_conformers().size(), 20 );
		TS_ASSERT_EQUALS( set6_->size(), 38 );
	}

	// Confirm that the proper conformers are loaded by get_ideal_conformer_by_name() and
	// get_ideal_conformer_by_CP_parameters.
	void test_get_ideal_conformer_by_descriptor_methods()
	{
		using namespace core;
		using namespace core::chemical::rings;
		using namespace utility;

		TR << "Testing get_ideal_conformer_by_name(), get_ideal_conformer_by_CP_parameters(), and "
			" get_ideal_conformer_from_nus() methods of RingConformerSet for 5- and 6-membered rings." << std::endl;

		// Set up variables.
		vector1< Real > params5, params6, params_bad;
		params5.resize( 2 );
		params6.resize( 3 );

		vector1< Angle > nus5, nus6, nus_bad;
		nus5.resize( 4 );
		nus6.resize( 5 );

		// Test some exact-match values.
		params5[ q ] = 0.4;
		params5[ PHI ] = 180.0;

		params6[ q ] = 0.55;
		params6[ PHI ] = 180.0;
		params6[ THETA ] = 90.0;

		nus5[ 1 ] = -30.0;
		nus5[ 2 ] = 0.0;
		nus5[ 3 ] = 30.0;
		nus5[ 4 ] = -60.0;

		nus6[ 1 ] = 0.0;
		nus6[ 2 ] = -60.0;
		nus6[ 3 ] = 60.0;
		nus6[ 4 ] = 0.0;
		nus6[ 5 ] = -60.0;

		TS_ASSERT_EQUALS( set5_->get_ideal_conformer_by_name( "EO" ).specific_name,
			set5_->get_ideal_conformer_by_CP_parameters( params5 ).specific_name );
		TS_ASSERT_EQUALS( set6_->get_ideal_conformer_by_name( "BO,3" ).specific_name,
			set6_->get_ideal_conformer_by_CP_parameters( params6 ).specific_name );

		TS_ASSERT_EQUALS( set5_->get_ideal_conformer_by_name( "EO" ).specific_name,
			set5_->get_ideal_conformer_from_nus( nus5 ).specific_name );
		TS_ASSERT_EQUALS( set6_->get_ideal_conformer_by_name( "BO,3" ).specific_name,
			set6_->get_ideal_conformer_from_nus( nus6 ).specific_name );

		// Test that rounding is handled properly.
		params5[ PHI ] = 43.2;  // should round to 36.0

		nus5[ 1 ] = 65.4;  // ideal is -60.0
		nus5[ 2 ] = -32.1;  // ideal is 30.0
		nus5[ 3 ] = -0.123;  // ideal is 0.0
		nus5[ 4 ] = 34.5;  // ideal is -30.0

		TS_ASSERT_EQUALS( set5_->get_ideal_conformer_by_name( "E1" ).specific_name,
			set5_->get_ideal_conformer_by_CP_parameters( params5 ).specific_name );
		TS_ASSERT_EQUALS( set5_->get_ideal_conformer_by_name( "E1" ).specific_name,
			set5_->get_ideal_conformer_from_nus( nus5 ).specific_name );

		params6[ PHI ] = 123.4;  // should round to 120.0
		params6[ THETA ] = 123.4;  // should round to 135.0

		nus6[ 1 ] = -0.987;  // ideal is 0.0
		nus6[ 2 ] = -0.123;  // ideal is 0.0
		nus6[ 3 ] = -32.1;  // ideal is -30.0
		nus6[ 4 ] = 67.8;  // ideal is 60.0
		nus6[ 5 ] = -67.8;  // ideal is -60.0

		TS_ASSERT_EQUALS( set6_->get_ideal_conformer_by_name( "5E" ).specific_name,
			set6_->get_ideal_conformer_by_CP_parameters( params6 ).specific_name );
		TS_ASSERT_EQUALS( set6_->get_ideal_conformer_by_name( "5E" ).specific_name,
			set6_->get_ideal_conformer_from_nus( nus6 ).specific_name );

		// Test non-principal angles.
		nus6[ 1 ] = -300.0;
		nus6[ 2 ] = 300.0;
		nus6[ 3 ] = 420.0;
		nus6[ 4 ] = -420.0;
		nus6[ 5 ] = 60.0;

		TS_ASSERT_EQUALS( set6_->get_ideal_conformer_by_name( "4C1" ).specific_name,
			set6_->get_ideal_conformer_from_nus( nus6 ).specific_name );

		// Test that chairs are handled correctly, as they reside at the poles where phi is meaningless.
		params6[ THETA ] = 0.0;

		TS_ASSERT_EQUALS( set6_->get_ideal_conformer_by_name( "4C1" ).specific_name,
			set6_->get_ideal_conformer_by_CP_parameters( params6 ).specific_name );

		// Test for bad input.
		try {
			set_throw_on_next_assertion_failure();
			set5_->get_ideal_conformer_by_name( "FOO" );  // This should force an exit.
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_Base const & e) {
			std::string expected_error_message( "ERROR: No conformer with given name found in this set; exiting.\n" );
			TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ), expected_error_message.size() ), expected_error_message );
		}

		try {
			set_throw_on_next_assertion_failure();
			set6_->get_ideal_conformer_by_CP_parameters( params_bad );  // This should force an exit.
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_Base const & e) {
			std::string expected( "ERROR: An N-membered ring is described by exactly N-3 Cremer-Pople parameters, "
				"yet a different number was provided; exiting.\n" );
			TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ), expected.size() ), expected );
		}
		params_bad.push_back( 0.0 );  // still bad because not enough params for a 6-membered ring
		try {
			set_throw_on_next_assertion_failure();
			set6_->get_ideal_conformer_by_CP_parameters( params_bad );  // This should force an exit.
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_Base const & e) {
			std::string expected( "ERROR: An N-membered ring is described by exactly N-3 Cremer-Pople parameters, "
				"yet a different number was provided; exiting.\n" );
			TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ), expected.size() ), expected );
		}
		params_bad.push_back( 180.0 );
		params_bad.push_back( 90.0 );  // still bad because q = 0.0 (planar)
		try {
			set_throw_on_next_assertion_failure();
			set6_->get_ideal_conformer_by_CP_parameters( params_bad );  // This should force an exit.
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_Base const & e) {
			std::string expected( "ERROR: Planar ring conformations are not handled by Rosetta; "
				"please specify a non-zero q value; exiting.\n" );
			TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ), expected.size() ), expected );
		}
		params_bad.push_back( 90.0 );  // still bad because too many params
		try {
			set_throw_on_next_assertion_failure();
			set6_->get_ideal_conformer_by_CP_parameters( params_bad );  // This should force an exit.
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_Base const & e) {
			std::string expected( "ERROR: An N-membered ring is described by exactly N-3 Cremer-Pople parameters, "
				"yet a different number was provided; exiting.\n" );
			TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ), expected.size() ), expected );
		}

		nus_bad.push_back( 180.0 );  // impossible angle; no such ring would be possible
		nus_bad.push_back( 180.0 );
		nus_bad.push_back( 180.0 );
		nus_bad.push_back( 180.0 );
		nus_bad.push_back( 180.0 );
		try {
			set_throw_on_next_assertion_failure();
			set6_->get_ideal_conformer_from_nus( nus_bad );  // This should force an exit.
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_Base const & e) {
			std::string expected( "ERROR: No conformer with given nu angles found in this set; exiting.\n" );
			TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ), expected.size() ), expected );
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

		TR << "Testing that RingConformer subsets have been populated correctly."  << std::endl;
		TS_ASSERT_EQUALS( set5_->get_lowest_energy_conformer().specific_name, "1E" );
		TS_ASSERT_EQUALS( set6_->get_lowest_energy_conformer().specific_name, "4C1" );

		vector1< string > expected_lows;
		expected_lows.push_back( "1E" );
		expected_lows.push_back( "2E" );
		expected_lows.push_back( "3E" );

		core::chemical::rings::RingConformer random_conformer_from_low_energy_subset;

		random_conformer_from_low_energy_subset = set5_->get_random_local_min_conformer();
		TS_ASSERT( expected_lows.contains( random_conformer_from_low_energy_subset.specific_name ) );

		expected_lows.push_back( "4C1" );
		random_conformer_from_low_energy_subset = set6_->get_random_local_min_conformer();
		TS_ASSERT( expected_lows.contains( random_conformer_from_low_energy_subset.specific_name ) );
	}

private:
	// Private data ////////////////////////////////////////////////////////////
	core::chemical::rings::RingConformerSetOP set5_;
	core::chemical::rings::RingConformerSetOP set6_;

};  // class RingConformerSetTests
