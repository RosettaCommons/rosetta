// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/trie/trie_vs_trie.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_etable_etrie_CountPairData_1_1_hh
#define INCLUDED_core_scoring_etable_etrie_CountPairData_1_1_hh

// Unit Headers
#include <core/scoring/etable/etrie/CountPairData_1_1.fwd.hh>

// Project Headers
#include <core/types.hh>

// STL Headers
#include <iosfwd>
#include <utility/assert.hh>

// Utility headers
#include <utility/assert.hh>

namespace core {
namespace scoring {
namespace etable {
namespace etrie {

class CountPairData_1_1
{
public:

	CountPairData_1_1() : connection_distance_( 0 ) {}

	void set_dist_to_connect_point( Size entry, Size connpoint, Size connection_dist );

	inline
	bool operator < ( CountPairData_1_1 const & other ) const
	{
		return connection_distance_ < other.connection_distance_;
	}

	inline
	bool operator == ( CountPairData_1_1 const & other ) const
	{
		return connection_distance_ == other.connection_distance_;
	}


	void print( std::ostream & os ) const;

	//void
	//set_count_pair_data_to_use( Size connection_id ) const;

	inline
	Size
	conn_dist( Size ASSERT_ONLY( which_connection) ) const
	{
	debug_assert( which_connection == 0 );
		return connection_distance_;
	}

private:
	Size connection_distance_;

};

std::ostream & operator << ( std::ostream & os, CountPairData_1_1 const & cpdat );

} // namespace etrie
} // namespace etable
} // namespace scoring
} // namespace core

#endif
