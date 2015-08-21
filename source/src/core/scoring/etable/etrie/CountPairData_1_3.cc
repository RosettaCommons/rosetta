// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/etrie/CountPairData_1_3.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/etable/etrie/CountPairData_1_3.hh>

// Utility headers
#include <utility/assert.hh>

// STL Headers
#include <iostream>
#include <utility/assert.hh>

namespace core {
namespace scoring {
namespace etable {
namespace etrie {

void CountPairData_1_3::set_dist_to_connect_point(
	Size entry,
	Size ASSERT_ONLY( connpoint ),
	Size connection_dist
)
{
	debug_assert( entry > 0 && entry <= 3 );
	debug_assert( connpoint == 1 );
	connection_distances_[ entry - 1 ] = connection_dist;
}


void CountPairData_1_3::print( std::ostream & os ) const
{
	os << "CountPairData_1_3" << std::endl;
}

/*
void
CountPairData_1_3::set_count_pair_data_to_use(
Size connection_id
) const
{
debug_assert( connection_id == 1 || connection_id == 2 || connection_id == 3 );
connection_distance_at_hand_ = connection_distances_[ connection_id - 1 ];
}
*/

std::ostream & operator << ( std::ostream & os, CountPairData_1_3 const & cpdat )
{
	cpdat.print( os );
	return os;
}


} // namespace etrie
} // namespace etable
} // namespace scoring
} // namespace core

