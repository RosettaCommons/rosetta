// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/etrie/CountPairData_1_3.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_etable_etrie_CountPairData_1_3_hh
#define INCLUDED_core_scoring_etable_etrie_CountPairData_1_3_hh

// Unit Headers
#include <core/scoring/etable/etrie/CountPairData_1_3.fwd.hh>

// Project Headers
#include <core/types.hh>

// STL Headers
#include <iosfwd>
#include <utility/assert.hh>

namespace core {
namespace scoring {
namespace etable {
namespace etrie {

class CountPairData_1_3
{
public:

	void set_dist_to_connect_point( Size entry, Size connpoint, Size connection_dist );

	void print( std::ostream & os ) const;

	inline
	bool operator < ( CountPairData_1_3 const & other ) const
	{
		if ( connection_distances_[ 0 ] < other.connection_distances_[ 0 ] ) return true;
		else if ( other.connection_distances_[ 0 ] < connection_distances_[ 0 ] ) return false;
		else if ( connection_distances_[ 1 ] < other.connection_distances_[ 1 ] ) return true;
		else if ( other.connection_distances_[ 1 ] < connection_distances_[ 1 ] ) return false;
		else if ( connection_distances_[ 2 ] < other.connection_distances_[ 2 ] ) return true;
		else if ( other.connection_distances_[ 2 ] < connection_distances_[ 2 ] ) return false;

		return false;
	}

	inline
	bool operator == ( CountPairData_1_3 const & other ) const
	{
		if ( connection_distances_[ 0 ] == other.connection_distances_[ 0 ] &&
				connection_distances_[ 1 ] == other.connection_distances_[ 1 ] &&
				connection_distances_[ 2 ] == other.connection_distances_[ 2 ] ) {
			return true;
		} else {
			return false;
		}
	}

	//void
	//set_count_pair_data_to_use( Size entry ) const;

	inline
	Size
	conn_dist( Size which_connection ) const
	{
		debug_assert( which_connection == 0 || which_connection == 1 || which_connection == 2 );
		return connection_distances_[ which_connection ];
	}

private:
	//mutable Size connection_distance_at_hand_;
	Size connection_distances_[ 3 ];

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

std::ostream & operator << ( std::ostream & os, CountPairData_1_3 const & cpdat );


} // namespace etrie
} // namespace etable
} // namespace scoring
} // namespace core

#endif
