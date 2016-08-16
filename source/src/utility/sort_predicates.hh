// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/
/// @brief  sort predicates for using std::pair in std::sort.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_utility_sort_predicates_hh
#define INCLUDED_utility_sort_predicates_hh

// STL headers
#include <utility> // std::pair

namespace utility {

template < class S, class T >
struct SortFirst
{
public:
	bool operator() ( std::pair< S, T > const & left, std::pair< S, T > const & right )
	{
		return left.first < right.first;
	}
};


template < class S, class T >
struct SortSecond
{
public:
	bool operator() ( std::pair< S, T > const & left, std::pair< S, T > const & right )
	{
		return left.second < right.second;
	}
};

}

#endif
