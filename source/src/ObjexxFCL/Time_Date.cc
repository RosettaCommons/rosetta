// Time and Date Functions
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// ObjexxFCL Headers
#include <ObjexxFCL/Time_Date.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/Fstring.hh>

// C++ Headers
#include <cassert>
#include <ctime>
#include <sstream>


namespace ObjexxFCL {


/// @brief Current Time: HH, MM, SS
void
itime( FArray1_int & timearray )
{
	assert( timearray.l() <= 1 );
	assert( timearray.u() >= 3 );
	std::time_t const current_time( std::time( nullptr ) );
	std::tm const * const timeinfo( std::localtime( &current_time ) );
	timearray( 1 ) = timeinfo->tm_hour;
	timearray( 2 ) = timeinfo->tm_min;
	timearray( 3 ) = timeinfo->tm_sec;
}


/// @brief Current Date: DD, MM, YYYY
void
idate( FArray1_int & datearray )
{
	assert( datearray.l() <= 1 );
	assert( datearray.u() >= 3 );
	std::time_t const current_time( std::time( nullptr ) );
	std::tm const * const timeinfo( std::localtime( &current_time ) );
	datearray( 1 ) = timeinfo->tm_mday; // Day of month: 1,2,...
	datearray( 2 ) = timeinfo->tm_mon + 1; // Month of year: 1-12
	datearray( 3 ) = timeinfo->tm_year + 1900; // Year
}


/// @brief Current Date Numeric (Not Y2K Compliant): MM, DD, YY
void
idate( int & month, int & day, int & year )
{
	std::time_t const current_time( std::time( nullptr ) );
	std::tm const * const timeinfo( std::localtime( &current_time ) );
	month = timeinfo->tm_mon + 1; // Month of year: 1-12
	day = timeinfo->tm_mday; // Day of month: 1,2,...
	year = timeinfo->tm_year - ( timeinfo->tm_year / 100 ) * 100; // Year: 0-99 (2-digit)
}


/// @brief Current Date String (Not Y2K Compliant): DD-MMM-YY
void
date( Fstring & day )
{
	assert( day.length() >= 9 );
	int m, d, y;
	idate( m, d, y );
	std::stringstream s;
	s << std::setfill( '0' ) << std::setw( 2 ) << d;
	day = "  -   -  ";
	day( 1, 2 ) = s.str();
	s.str( "" );
	s << std::setfill( '0' ) << std::setw( 2 ) << y;
	day( 8, 9 ) = s.str();
	Fstring mmm( 3 );
	switch ( m ) {
	case 1 :
		mmm = "JAN";
		break;
	case 2 :
		mmm = "FEB";
		break;
	case 3 :
		mmm = "MAR";
		break;
	case 4 :
		mmm = "APR";
		break;
	case 5 :
		mmm = "MAY";
		break;
	case 6 :
		mmm = "JUN";
		break;
	case 7 :
		mmm = "JUL";
		break;
	case 8 :
		mmm = "AUG";
		break;
	case 9 :
		mmm = "SEP";
		break;
	case 10 :
		mmm = "OCT";
		break;
	case 11 :
		mmm = "NOV";
		break;
	case 12 :
		mmm = "DEC";
		break;
	default :
		assert( false );
		break;
	}
	day( 4, 6 ) = mmm;
}


} // namespace ObjexxFCL
