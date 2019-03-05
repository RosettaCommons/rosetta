// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    utility/dating.cc
/// @brief   Definitions for utility functions involving calendar dates, not for
///          finding new romantic partners.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @note    This file was created at the request of Vikram.  Currently, the only
///          code using this (currently) single function is .pdb file output
///          code, but perhaps other needs for dating output might be needed in
///          the future....


// Unit Headers
#include <utility/dating.hh>

// C++ Headers
#include <sstream>


namespace utility {

/// @brief  Return current date in the requested format.
std::string
get_current_date( DateFormat const format )
{
	using namespace std;

	time_t const t( time( nullptr ) );
	tm date_and_time( *localtime( &t ) );
	stringstream day;
	string month;
	stringstream year;
	string date;

	if ( format == PDB_FORMAT ) {
		// These formats use a leading zero.
		day << ( date_and_time.tm_mday < 10 ? "0" : "" ) << date_and_time.tm_mday;
	}

	if ( format == PDB_FORMAT ) {
		// These formats use 3-letter, all-caps month names.
		switch ( date_and_time.tm_mon ) {
		case 0 :
			month = "JAN";
			break;
		case 1 :
			month = "FEB";
			break;
		case 2 :
			month = "MAR";
			break;
		case 3 :
			month = "APR";
			break;
		case 4 :
			month = "MAY";
			break;
		case 5 :
			month = "JUN";
			break;
		case 6 :
			month = "JUL";
			break;
		case 7 :
			month = "AUG";
			break;
		case 8 :
			month = "SEP";
			break;
		case 9 :
			month = "OCT";
			break;
		case 10 :
			month = "NOV";
			break;
		default :  // case 11
			month = "DEC";
		}
	}

	// System dates are in years since 1900.
	if ( date_and_time.tm_year > 99 ) {
		date_and_time.tm_year -= 100;
	}
	if ( format == PDB_FORMAT ) {
		// These formats use 2-letter years.
		year << ( date_and_time.tm_year < 10 ? "0" : "" ) << date_and_time.tm_year;
	}

	if ( format == PDB_FORMAT ) {
		// dd-MM-yy
		date = day.str() + "-" + month + "-" + year.str();
	}

	return date;
}


}  // namespace utility
