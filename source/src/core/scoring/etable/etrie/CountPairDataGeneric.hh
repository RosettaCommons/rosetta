// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/etrie/CountPairDataGeneric.hh
/// @brief  Class to hold per-atom count pair for residues with non-canonical inter-residue connections
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_etable_etrie_CountPairDataGeneric_hh
#define INCLUDED_core_scoring_etable_etrie_CountPairDataGeneric_hh

// Unit Headers
#include <core/scoring/etable/etrie/CountPairDataGeneric.fwd.hh>

// Project Headers
#include <core/types.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector0.hh>
#include <utility/vector1.hh>

// STL Headers
#include <iosfwd>

#include <platform/types.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <utility/assert.hh>
#include <cstddef>
#include <limits>
#include <vector>


namespace core {
namespace scoring {
namespace etable {
namespace etrie {

class GenericResidueConnectionData
{
private:
   utility::vector1< Size > path_distances_to_connection_points_;

public:

	GenericResidueConnectionData();

   GenericResidueConnectionData( Size n_connections );

	inline
	bool operator < ( GenericResidueConnectionData const & other ) const
	{
	debug_assert( size() == other.size() );
		for ( Size ii = 1; ii <= path_distances_to_connection_points_.size(); ++ii ) {
			if ( path_distances_to_connection_points_[ ii ] < other.path_distances_to_connection_points_[ ii ] ) {
				return true;
			} else if ( other.path_distances_to_connection_points_[ ii ] != path_distances_to_connection_points_[ ii ] ) {
				return false;
			}
		}
		return false;
	}

	inline
	bool operator == ( GenericResidueConnectionData const & other ) const
	{
	debug_assert( size() == other.size() );
		for ( Size ii = 1; ii <= path_distances_to_connection_points_.size(); ++ii ) {
			if ( path_distances_to_connection_points_[ ii ] != other.path_distances_to_connection_points_[ ii ] ) {
				return false;
			}
		}
		return true;
	}

	inline
	bool operator != ( GenericResidueConnectionData const & other ) const
	{
		return ! operator == ( other );
	}


	/// @brief getter
	Size
	operator [] ( Size connection ) const {
		return path_distances_to_connection_points_[ connection ];
	}

	void set_dist_to_connect_point( Size connpoint, Size connection_dist );

	Size
	size() const {
		return path_distances_to_connection_points_.size();
	}
};

class CountPairDataGeneric
{
private:

	/// Vector0 to match index-from-0 convention that the other CountPairData classes use (since they hold c-style arrays).
	///iwd  Lies!  All other CountPair-related classes use indexing from 1.  The above comment appears to be out-of-date.
	/// Almost all code for this class assumes indexing from 1.
	utility::vector1< GenericResidueConnectionData > residue_connection_data_;
	//mutable GenericResidueConnectionData const * data_at_hand_;

public:
	CountPairDataGeneric();

	/// @brief "entry" is indexed from 1
	void set_dist_to_connect_point( Size entry, Size connpoint, Size connection_dist );

	void print( std::ostream & os ) const;

	inline
	bool operator < ( CountPairDataGeneric const & other ) const
	{
		for( Size ii = 1; ii <= residue_connection_data_.size(); ++ii ) {
			if ( residue_connection_data_[ ii ] < other.residue_connection_data_[ ii ] ) {
				return true;
			} else if ( other.residue_connection_data_[ ii ] != residue_connection_data_[ ii ] ) {
				return false;
			}
		}
		return false;
	}

	inline
	bool operator == ( CountPairDataGeneric const & other ) const
	{

		for (Size ii = 1; ii <= residue_connection_data_.size(); ++ii ) {
			if (residue_connection_data_[ ii ] != other.residue_connection_data_[ ii ])
				return false;
		}
		return true;
	}

	//void
	//set_count_pair_data_to_use( Size entry ) const;

	/// For use by the 1-connection TrieCountPair functions
	/// Indexed from 0.
	inline
	Size
	conn_dist( Size which_connection ) const
	{
	debug_assert ( residue_connection_data_[ which_connection+1 ].size() == 1 );
		return residue_connection_data_[ which_connection+1 ][1];
	}

	inline
	Size
	n_connected_residues() const {
		return residue_connection_data_.size();
	}

	/// Indexed from 0.
	inline
	GenericResidueConnectionData const &
	conn_dat( Size which_connection ) const {
		return residue_connection_data_[ which_connection+1 ];
	}

};

std::ostream & operator << ( std::ostream & os, CountPairDataGeneric const & cpdat );


} // namespace etrie
} // namespace etable
} // namespace scoring
} // namespace core

#endif
