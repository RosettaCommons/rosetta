// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
// This file is based off of code from the Biochemistry Library (BCL).
// The BCL is copyright Vanderbilt University (Meiler Lab), a RosettaCommons member

/// @file   core/chemical/gasteiger/util.hh
/// @brief  Utilities for dealing with gasteiger things.
/// @author To Rosetta transitioning: Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_gasteiger_util_hh
#define INCLUDED_core_chemical_gasteiger_util_hh

#include <numeric/util.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>

#include <set>
#include <iostream>

namespace core {
namespace chemical {
namespace gasteiger {

inline
void safe_write( std::ostream & out, core::Size const &val, bool sep = true ) {
	if ( numeric::is_undefined(val) ) {
		out << "nan";
	} else {
		out << val;
	}
	if ( sep ) {
		out << ' ';
	}
}

inline
void safe_write( std::ostream & out, core::Real const &val, bool sep = true ) {
	if ( numeric::is_undefined(val) ) {
		out << "nan";
	} else {
		out << val;
	}
	if ( sep ) {
		out << ' ';
	}
}

inline
void safe_read( std::istream & in, core::Size & val ) {
	in >> val;
	if( ! in.good() ) { // Error with input (e.g. not an int)
		in.clear();  // Reset input error flags
		std::string tag;
		in >> tag;
		if ( tag == "nan" ) {
			val = numeric::get_undefined_size();
		} else {
			utility_exit_with_message( "Malformatted file - got alphabetic when expecting numeric.");
		}
	}
}

inline
void safe_read( std::istream & in, core::Real & val ) {
	in >> val;
	if( ! in.good() ) { // Error with input (e.g. not an int)
		in.clear();  // Reset input error flags
		std::string tag;
		in >> tag;
		if ( tag == "nan" ) {
			val = numeric::get_undefined_real();
		} else {
			utility_exit_with_message( "Malformatted file - got alphabetic when expecting numeric.");
		}
	}
}

template < typename T >
std::set< T > parse_enum_set( std::string in ) {
	//The first two letters are an abbreviation of which enum it is -- this should be pre-checked.
	std::set< T > out;
	for( core::Size ii(2); ii < in.size(); ++ii ) {
		char c( in[ii] );
		int val;
		// Assuming reasonable code page layout where 0-9 and A-Z are each consecutive ascending in their range
		if ( c >= '0' && c <= '9' ) {
			val = c - '0';
		} else if ( c >= 'A' && c <= 'Z' ) {
			val = c - 'A' + 10;
		} else {
			utility_exit_with_message( "Unknown compact enum value: "+c );
		}
		out.insert( T( val ) );
	}
	return out;
}

template < typename T >
std::string compact_enum_set( std::set< T > in, std::string prefix = "xx" ) {
	//The first two letters are an abbreviation of which enum it is -- this should be pre-checked.
	std::string out(prefix);
	for( typename std::set< T >::iterator iter( in.begin() ), end( in.end() ); iter != end; ++iter) {
		int val( *iter );
		// Assuming reasonable code page layout where 0-9 and A-Z are each consecutive ascending in their range
		if ( val >= 0 && val <= 9 ) {
			out.push_back( '0' + val );
		} else if ( val >= 10 && val <= 36 ) {
			out.push_back( 'A' + val - 10 );
		} else {
			utility_exit_with_message( "Encountered uncompactable enum value." );
		}
	}
	return out;
}




} // namespace gasteiger
} // namespace core
} // namespace chemical

  #endif
