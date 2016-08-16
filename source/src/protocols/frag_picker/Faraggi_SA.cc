// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/Faraggi_SA.cc
/// @brief  Faraggi maximum SA values from Faraggi et al. Proteins 2008
/// @author David E Kim

#include <protocols/frag_picker/Faraggi_SA.hh>

// utility headers
#include <utility/exit.hh>

// project headers

// C++ headers
#include <map>

namespace protocols {
namespace frag_picker {

// BEGIN local functions

/// @brief setup the faraggi SA max map, see Faraggi et al. Proteins 2008
std::map< char, core::Real > setup_faraggi_map() {
	std::map< char, core::Real > sa_max;
	sa_max[ 'R' ] = 271.0;
	sa_max[ 'K' ] = 257.0;
	sa_max[ 'D' ] = 183.0;
	sa_max[ 'E' ] = 286.0;
	sa_max[ 'N' ] = 188.0;
	sa_max[ 'Q' ] = 215.0;
	sa_max[ 'H' ] = 238.0;
	sa_max[ 'Y' ] = 250.0;
	sa_max[ 'W' ] = 260.0;
	sa_max[ 'S' ] = 181.0;
	sa_max[ 'T' ] = 192.0;
	sa_max[ 'G' ] = 136.0;
	sa_max[ 'P' ] = 170.0;
	sa_max[ 'A' ] = 169.0;
	sa_max[ 'M' ] = 236.0;
	sa_max[ 'C' ] = 139.0;
	sa_max[ 'F' ] = 221.0;
	sa_max[ 'L' ] = 221.0;
	sa_max[ 'V' ] = 171.0;
	sa_max[ 'I' ] = 210.0;
	sa_max[ 'r' ] = 271.0;
	sa_max[ 'k' ] = 257.0;
	sa_max[ 'd' ] = 183.0;
	sa_max[ 'e' ] = 286.0;
	sa_max[ 'n' ] = 188.0;
	sa_max[ 'q' ] = 215.0;
	sa_max[ 'h' ] = 238.0;
	sa_max[ 'y' ] = 250.0;
	sa_max[ 'w' ] = 260.0;
	sa_max[ 's' ] = 181.0;
	sa_max[ 't' ] = 192.0;
	sa_max[ 'g' ] = 136.0;
	sa_max[ 'p' ] = 170.0;
	sa_max[ 'a' ] = 169.0;
	sa_max[ 'm' ] = 236.0;
	sa_max[ 'c' ] = 139.0;
	sa_max[ 'f' ] = 221.0;
	sa_max[ 'l' ] = 221.0;
	sa_max[ 'v' ] = 171.0;
	sa_max[ 'i' ] = 210.0;
	return sa_max;
}

/// @brief faraggi SA max map
inline std::map< char, core::Real > & sa_faraggi_max() {
	// static initialization only happens once
	static std::map< char, core::Real > * sa_faraggi_max_ = new std::map< char, core::Real >( setup_faraggi_map() );
	return *sa_faraggi_max_;
}

// END local functions

core::Real
sa_faraggi_max( char aa )
{
	std::map< char, core::Real >::const_iterator iter = sa_faraggi_max().find( aa );
	if ( iter == sa_faraggi_max().end() ) {
		utility_exit_with_message( "unrecognized sa_faraggi_max aa " + std::string(1, aa) );
	}
	return iter->second;
}

} // frag_picker
} // protocols


