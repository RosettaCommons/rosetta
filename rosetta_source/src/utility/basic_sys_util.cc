// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/basic_sys_util.cc
/// @brief  All system functions in utility that have no other home
/// @author David Kim (dekim@u.washington.edu)
/// @todo   Break out platform-specific code.
/// @todo   Get rid of output messages: unnecessary dependence on <iostream>


// Unit headers
#include <utility/basic_sys_util.hh>

// C++ headers
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

// Platform headers
#ifdef _WIN32
#include <windows.h>
#include <float.h>
#else // Not _WIN32
#include <unistd.h>
#endif // _WIN32

// Boinc headers
#ifdef BOINC
#include <utility/boinc/boinc_util.hh>
#endif


namespace utility {


/// @brief Sleep for a specified number of seconds
void
sys_sleep( double const seconds )
{
	using std::cout;
	using std::endl;
	using std::fmod;

	//cout << "sleep" << endl;

#ifdef _WIN32
	::Sleep( (int) (1000 * seconds) );
#else // Not _WIN32
	unsigned int remaining_time = (int) seconds;
	while ( true ) {
		remaining_time = ::sleep( remaining_time );
		if ( remaining_time == 0 ) break;
		if ( remaining_time > seconds ) break; // paranoia
	}
	int const x = static_cast< int >( fmod( seconds * 1000000, 1000000 ) );
	if ( x ) ::usleep(x);
#endif // _WIN32

	//cout << "sleep end" << endl;
}


/// @brief Generate timestamp string
std::string
timestamp()
{
	using std::ostringstream;
	using std::setw;
	using std::time;
	using std::time_t;
	using std::tm;

	time_t currentTime = time( 0 );
	struct tm * now = std::localtime( &currentTime );

	ostringstream timestamp;
	timestamp
	 << "["
	 << setw( 4 ) << ( now->tm_year + 1900 ) << "-"
	 << setw( 2 ) << ( now->tm_mon + 1 ) << "-"
	 << setw( 2 ) << ( now->tm_mday ) << " "
	 << setw( 2 ) << ( now->tm_hour ) << ":"
	 << setw( 2 ) << ( now->tm_min ) << ":"
	 << setw( 2 ) << ( now->tm_sec ) << ":"
	 << "]";

	return timestamp.str();
}


} // namespace utility
