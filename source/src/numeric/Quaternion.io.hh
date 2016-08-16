// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/Quaternion.io.hh
/// @brief  Quaternion input/output functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_Quaternion_io_hh
#define INCLUDED_numeric_Quaternion_io_hh


// Package headers
#include <numeric/Quaternion.hh>
#include <numeric/IOTraits.hh>

// C++ headers
#include <iostream>
#include <iomanip>
#include <sstream>


namespace numeric {


/// @brief stream << Quaternion output operator
template< typename T >
std::ostream &
operator <<( std::ostream & stream, Quaternion< T > const & q )
{
	// Types
	using std::setw;
	typedef  IOTraits< T >  Traits;

	// Save current stream state and set persistent state
	std::ios_base::fmtflags const old_flags = stream.flags();
	int const old_precision = stream.precision( Traits::precision() );
	stream << std::right << std::showpoint << std::uppercase;

	// Output Quaternion
	int const w = Traits::width();
	stream << setw( w ) << q.w() << ' ' << setw( w ) << q.x() << ' ' << setw( w ) << q.y() << ' ' << setw( w ) << q.z();

	// Restore previous stream state
	stream.precision( old_precision );
	stream.flags( old_flags );

	return stream;
}


/// @brief stream >> Quaternion input operator
/// @note Supports whitespace-separated values with optional commas between values as
///       long as whitespace is also present
/// @note Quaternion can optionally be enclosed in parentheses () or square brackets []
/// @note String or char values containing whitespace or commas or enclosed in quotes
///       are not supported
template< typename T >
std::istream &
operator >>( std::istream & stream, Quaternion< T > & q )
{
	bool parens = false; // Opening ( present?
	bool brackets = false; // Opening [ present?

	{ // w
		std::string input_string;
		stream >> input_string;
		if ( input_string == "(" ) { // Skip opening (
			stream >> input_string;
			parens = true;
		} else if ( input_string[ 0 ] == '(' ) { // Skip opening (
			input_string.erase( 0, 1 );
			brackets = true;
		} else if ( input_string == "[" ) { // Skip opening [
			stream >> input_string;
			brackets = true;
		} else if ( input_string[ 0 ] == '[' ) { // Skip opening [
			input_string.erase( 0, 1 );
			brackets = true;
		}
		std::string::size_type const input_size = input_string.size();
		if ( ( input_size > 0 ) && ( input_string[ input_size - 1 ] == ',' ) ) {
			input_string.erase( input_size - 1 ); // Remove trailing ,
		}
		std::istringstream num_stream( input_string );
		num_stream >> q.w();
	}

	{ // x
		std::string input_string;
		stream >> input_string;
		if ( input_string == "," ) { // Skip ,
			stream >> input_string;
		} else if ( input_string[ 0 ] == ',' ) { // Skip leading ,
			input_string.erase( 0, 1 );
		}
		std::string::size_type const input_size = input_string.size();
		if ( ( input_size > 0 ) && ( input_string[ input_size - 1 ] == ',' ) ) {
			input_string.erase( input_size - 1 ); // Remove trailing ,
		}
		std::istringstream num_stream( input_string );
		num_stream >> q.x();
	}

	{ // y
		std::string input_string;
		stream >> input_string;
		if ( input_string == "," ) { // Skip ,
			stream >> input_string;
		} else if ( input_string[ 0 ] == ',' ) { // Skip leading ,
			input_string.erase( 0, 1 );
		}
		std::string::size_type const input_size = input_string.size();
		if ( ( input_size > 0 ) && ( input_string[ input_size - 1 ] == ',' ) ) {
			input_string.erase( input_size - 1 ); // Remove trailing ,
		}
		std::istringstream num_stream( input_string );
		num_stream >> q.y();
	}

	{ // z
		std::string input_string;
		stream >> input_string;
		if ( input_string == "," ) { // Skip ,
			stream >> input_string;
		} else if ( input_string[ 0 ] == ',' ) { // Skip leading ,
			input_string.erase( 0, 1 );
		}
		std::string::size_type input_size = input_string.size();
		if ( parens || brackets ) { // Remove closing ) or ]
			if ( input_size > 0 ) {
				if ( parens ) {
					if ( input_string[ input_size - 1 ] == ')' ) { // Remove closing )
						input_string.erase( input_size - 1 );
						--input_size;
					}
				} else if ( brackets ) {
					if ( input_string[ input_size - 1 ] == ']' ) { // Remove closing ]
						input_string.erase( input_size - 1 );
						--input_size;
					}
				}
			}
		}
		if ( ( input_size > 0 ) && ( input_string[ input_size - 1 ] == ',' ) ) {
			input_string.erase( input_size - 1 ); // Remove trailing ,
		}
		std::istringstream num_stream( input_string );
		num_stream >> q.z();
	}

	// Remove closing ) or ] if opening ( or [ present
	if ( parens || brackets ) { // Remove closing ) or ]
		while ( ( stream.peek() == ' ' ) || ( stream.peek() == '\t' ) ) {
			stream.ignore();
		}
		if ( parens ) { // Remove closing ) if present
			if ( stream.peek() == ')' ) stream.ignore();
		} else if ( brackets ) { // Remove closing ] if present
			if ( stream.peek() == ']' ) stream.ignore();
		}
	}

	return stream;
}


} // namespace numeric


#endif // INCLUDED_numeric_Quaternion_io_HH
