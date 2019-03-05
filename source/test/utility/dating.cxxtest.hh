// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/utility/dating.cxxtest.hh
/// @brief   Test suite for utility functions involving calendar dates, not for
///          finding new romantic partners.
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
//#include <test/core/init_util.hh>

// Unit header
#include <utility/dating.hh>

// Basic header
#include <basic/Tracer.hh>

// C++ header
#include <string>


class DatingUtilitiesTests : public CxxTest::TestSuite {
public: // Standard methods ////////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		//core_init();
	}

	// Destruction
	void tearDown()
	{}


public: // Tests //////////////////////////////////////////////////////////////
	// Confirm that the returned date makes sense.
	void test_get_current_date()
	{
		using namespace std;
		using namespace utility;

		bool day_is_valid( false ), month_is_valid( false ), year_is_valid( false );
		string const date( get_current_date( PDB_FORMAT ) );
		int const day( stoi( date.substr( 0, 2 ) ) );
		string const month( date.substr( 3, 3 ) );
		int const year( stoi( date.substr( 7, 2 ) ) );

		if ( ( year > 0 ) && ( year <= 99 ) ) {
			year_is_valid = true;
		}

		if ( ( month == "JAN" ) ||
				( month == "FEB" ) ||
				( month == "MAR" ) ||
				( month == "APR" ) ||
				( month == "MAY" ) ||
				( month == "JUN" ) ||
				( month == "JUL" ) ||
				( month == "AUG" ) ||
				( month == "SEP" ) ||
				( month == "OCT" ) ||
				( month == "NOV" ) ||
				( month == "DEC" ) ) {
			month_is_valid = true;
		}

		if ( ( day > 0 ) ) {
			if ( month == "FEB" ) {
				if ( day <= 29 ) { day_is_valid = true; }
			} else if ( ( month == "SEP" ) || ( month == "APR" ) || ( month == "JUN" ) || ( month == "NOV" ) ) {
				if ( day <= 30 ) { day_is_valid = true; }
			} else {
				if ( day <= 31 ) { day_is_valid = true; }
			}
		}

		TS_ASSERT( day_is_valid );
		TS_ASSERT( month_is_valid );
		TS_ASSERT( year_is_valid );
	}
};
