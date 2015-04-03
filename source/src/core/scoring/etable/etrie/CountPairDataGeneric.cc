// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/etrie/CountPairDataGeneric.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/etable/etrie/CountPairDataGeneric.hh>

// Utility headers

// STL Headers
#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace etable {
namespace etrie {

GenericResidueConnectionData::GenericResidueConnectionData() {}

GenericResidueConnectionData::GenericResidueConnectionData( Size n_connections ) :
	path_distances_to_connection_points_( n_connections, 0 )
{}


void
GenericResidueConnectionData::set_dist_to_connect_point(
	Size connpoint,
	Size connection_dist
)
{
debug_assert( connpoint > 0 );
	if ( connpoint > path_distances_to_connection_points_.size()) {
		path_distances_to_connection_points_.resize( connpoint );
	}
	path_distances_to_connection_points_[ connpoint ] = connection_dist;
}


void CountPairDataGeneric::set_dist_to_connect_point(
	Size entry,
	Size connpoint,
	Size connection_dist
)
{
debug_assert( entry > 0);
	if ( entry > residue_connection_data_.size() ) {
		residue_connection_data_.resize(entry);
	}
	residue_connection_data_[ entry ].set_dist_to_connect_point( connpoint, connection_dist );
}

CountPairDataGeneric::CountPairDataGeneric(): residue_connection_data_(4)
{
	//std::cout << "Constructed  Count Pair Data Generic " << std::endl;
}

void CountPairDataGeneric::print( std::ostream & os ) const
{
	os << "CountPairDataGeneric, " << residue_connection_data_.size() << " conections:\n";
	for ( Size ii = 1; ii <= residue_connection_data_.size(); ++ii ) {
	   os << ii << " :";
		for ( Size jj = 1; jj <= residue_connection_data_[ ii ].size(); ++jj ) {
			os << " " << residue_connection_data_[ ii ][ jj ];
		}
		os << "\n";
	}
}


/*void
CountPairDataGeneric::set_count_pair_data_to_use(
	Size connection_id
) const
{
debug_assert( connection_id > 0 && connection_id <= residue_connection_data_.size());
	data_at_hand_ = &residue_connection_data_[ connection_id  ];
}
*/

std::ostream & operator << ( std::ostream & os, CountPairDataGeneric const & cpdat )
{
	cpdat.print( os );
	return os;
}


} // namespace etrie
} // namespace etable
} // namespace scoring
} // namespace core

