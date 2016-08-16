// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file utility/map_util.hh
/// @brief Utility functions for std::maps
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_utility_util_hh
#define INCLUDED_utility_util_hh

#include <map>

namespace utility {


/// @brief Does the map have the key?
template < class T >
bool
has_key(std::map< T, T > const & a_map, T const & key);


} //utility


#endif	//utility_util_hh

