// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/kdtree/KDPointList.cc
/// @brief definition for a class that keeps track of the best N KDPoint objects
/// by distance.
/// @author James Thompson

#include <numeric/types.hh>

#include <numeric/kdtree/KDPoint.hh>
#include <numeric/kdtree/KDPointList.hh>
#include <numeric/kdtree/constants.hh>

#include <utility/vector1.hh>

#include <iostream>
#include <algorithm>

namespace numeric {
namespace kdtree {

/// @details Auto-generated virtual destructor
KDPointList::~KDPointList() = default;

bool KDPoint_MinDist::operator() ( KDPointOP const & a, KDPointOP const & b ) {
	return ( a->distance() < b->distance() );
}

utility::vector1< KDPointOP >::const_iterator KDPointList::begin() const {
	return container_.begin();
}

utility::vector1< KDPointOP >::const_iterator KDPointList::end() const {
	return container_.end();
}

utility::vector1< KDPointOP >::iterator KDPointList::begin() {
	return container_.begin();
}

utility::vector1< KDPointOP >::iterator KDPointList::end() {
	return container_.end();
}

KDPointList::KDPointList( numeric::Size const n_to_keep ) :
	max_vals_( n_to_keep ),
	distance_cutoff_( REALLY_BIG_DISTANCE )
{}

void KDPointList::insert( KDPointOP pt ) {
	update_heap_();
	container_.push_back( pt );
	std::push_heap( container_.begin(), container_.end(), KDPoint_MinDist() );
	update_size_();
}

numeric::Real KDPointList::worst_distance() const {
	Real worst_dist = distance_cutoff();
	if ( container_.size() == max_values() ) {
		worst_dist = std::min( container_.front()->distance(), worst_dist );
	}

	return worst_dist;
}

numeric::Real KDPointList::distance_cutoff() const {
	return distance_cutoff_;
}

void KDPointList::distance_cutoff( numeric::Real const cutoff ) {
	distance_cutoff_ = cutoff;
}

KDPointOP KDPointList::operator[]( numeric::Size const pos ) const {
	assert( pos > 0 );
	assert( pos <= size() );
	return container_[pos];
}

numeric::Size KDPointList::size() const {
	return container_.size();
}

numeric::Size KDPointList::max_values() const {
	return max_vals_;
}

utility::vector1< KDPointOP > KDPointList::sorted_values() {
	using utility::vector1;
	vector1< KDPointOP > values( container_.begin(), container_.end() );
	sort_heap( values.begin(), values.end(), KDPoint_MinDist() );
	return values;
}

/// @brief merge another KDPointList with this KDPointList.
/// This calls insert which is a little slow, and is a candidate for
/// optimization if the insert() method shows up in profiling.
void KDPointList::merge( KDPointList const & other ) {
	for ( auto tmp : other ) {
		insert( tmp );
	}
}

void KDPointList::show( std::ostream & out ) const {
	out << "KDPointList with " << size() << " points." << std::endl;
	 for ( auto const & it : container_ ) {
		it->show( out );
		out << std::endl;
	}
}

// private function declarations
void KDPointList::update_heap_() {
	make_heap( container_.begin(), container_.end(), KDPoint_MinDist() );
}

void KDPointList::update_size_() {
	while ( container_.size() > max_vals_ ) {
		std::pop_heap( container_.begin(), container_.end(), KDPoint_MinDist() );
		container_.pop_back();
	}
}

} // kdtree
} // numeric
