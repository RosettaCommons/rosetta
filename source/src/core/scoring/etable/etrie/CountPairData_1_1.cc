// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/etrie/CountPairData_1_1.cc
/// @brief  CountPair Data for a residue with one connection point
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/etable/etrie/CountPairData_1_1.hh>

#include <utility/assert.hh>

// STL Headers
#include <utility/assert.hh>
#include <iostream>

namespace core {
namespace scoring {
namespace etable {
namespace etrie {

void CountPairData_1_1::set_dist_to_connect_point(
	Size ASSERT_ONLY(entry),
	Size ASSERT_ONLY(connpoint),
	Size connection_dist
)
{
	debug_assert( entry == 1 );
	debug_assert( connpoint == 1 );
	connection_distance_ = connection_dist;
}

void
CountPairData_1_1::print( std::ostream & os ) const
{
	os << "CountPairData_1_1: " << connection_distance_ << std::endl;
}

/*
void
CountPairData_1_1::set_count_pair_data_to_use(
Size ASSERT_ONLY(entry)
) const
{
debug_assert( entry == 0 );
}
*/

std::ostream & operator << ( std::ostream & os, CountPairData_1_1 const & cpdat )
{
	cpdat.print( os );
	return os;
}

} // namespace etrie
} // namespace etable
} // namespace scoring
} // namespace core

