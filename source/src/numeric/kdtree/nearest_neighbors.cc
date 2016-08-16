// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/kdtree/util.cc
/// @brief
/// @author James Thompson

#include <numeric/types.hh>

#include <numeric/kdtree/constants.hh>
#include <numeric/kdtree/nearest_neighbors.hh>
#include <numeric/kdtree/util.hh>
#include <numeric/kdtree/calc_distance.hh>
#include <numeric/kdtree/KDNode.hh>
#include <numeric/kdtree/KDTree.hh>
#include <numeric/kdtree/HyperRectangle.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>

#include <cmath>
#include <algorithm>

namespace numeric {
namespace kdtree {

void
nearest_neighbor(
	KDTree & tree,
	utility::vector1< numeric::Real > const & pt,
	// returns:
	KDNodeOP & nearest,
	numeric::Real & dist_sq
) {
	HyperRectangle bounds = (*tree.bounds());
	//numeric::Real max_dist_sq( sq_vec_distance(
	// bounds->lower(),
	// bounds->upper()
	//) );
	numeric::Real max_dist_sq( REALLY_BIG_DISTANCE );
	KDNodeOP root( tree.root() );

	nearest_neighbor(
		root, pt, bounds, max_dist_sq, nearest, dist_sq
	);
}

KDPointList
nearest_neighbors(
	KDTree & tree,
	utility::vector1< numeric::Real > const & pt,
	Size const wanted
) {
	HyperRectangle bounds = (*tree.bounds());
	//numeric::Real max_dist_sq( sq_vec_distance(
	// bounds->lower(),
	// bounds->upper()
	//) );
	numeric::Real max_dist_sq( REALLY_BIG_DISTANCE );
	KDNodeOP root( tree.root() );

	KDPointList nearest( wanted );

	nearest_neighbors(
		root, pt, bounds, max_dist_sq, nearest
	);

	return nearest;
}

KDPointList
nearest_neighbors(
	KDTree & tree,
	utility::vector1< numeric::Real > const & pt,
	Size const wanted,
	numeric::Real const max_dist_allowed
) {
	HyperRectangle bounds = (*tree.bounds());
	numeric::Real max_dist_sq( max_dist_allowed * max_dist_allowed );
	KDNodeOP root( tree.root() );

	KDPointList nearest( wanted );
	nearest.distance_cutoff( max_dist_sq );

	//std::cout << "distance cutoff = " << max_dist_sq << std::endl;

	nearest_neighbors(
		root, pt, bounds, max_dist_sq, nearest
	);

	return nearest;
}

void
nearest_neighbors(
	KDNodeOP & current,
	utility::vector1< numeric::Real > const & pt,
	HyperRectangle & bounds,
	numeric::Real max_dist_sq,

	// returns:
	KDPointList & neighbors
) {
	using numeric::Size;
	using numeric::Real;
	using utility::vector1;

	if ( !current ) {
		return;
	}
	Size const split_axis ( current->split_axis() );
	Real const split_value( current->location()[split_axis] );

	//std::cout << "currently at tree " << std::endl;
	//print_tree( std::cout, current, 1 );
	//std::cout << "max_dist_sq = " << max_dist_sq << std::endl;
	//std::cout << "split on ";
	//print_point( std::cout, pt );
	//std::cout << " (axis,value) = " << split_axis << "," << split_value << std::endl;

	// mid is a hyper-plane through this location and perpendicular
	// to the split axis
	vector1< Real >
		left_lower ( bounds.lower() ), left_upper ( bounds.upper() ),
		right_lower( bounds.lower() ), right_upper( bounds.upper() );
	left_upper [split_axis] = split_value;
	right_lower[split_axis] = split_value;

	//std::cout << "left ";
	HyperRectangle left_hr ( left_upper,  left_lower );
	//std::cout << "right ";
	HyperRectangle right_hr( right_upper, right_lower );

	KDNodeOP nearer, further;
	HyperRectangle nearer_hr, further_hr;

	if ( pt[split_axis] <= split_value ) {
		//std::cout << "chose left as nearer" << std::endl;
		nearer     = current->left_child();
		further    = current->right_child();
		nearer_hr  = left_hr;
		further_hr = right_hr;
	} else {
		//std::cout << "chose right as nearer" << std::endl;
		nearer     = current->right_child();
		further    = current->left_child();
		nearer_hr  = left_hr;
		further_hr = right_hr;
	}

	//std::cout << "calling nearest_neighbors on nearer child" << std::endl;

	nearest_neighbors(
		nearer, pt, nearer_hr, max_dist_sq, neighbors
	);
	max_dist_sq = std::min( neighbors.worst_distance(), max_dist_sq );

	// we need to seach this point and right-hand child (furthest
	// hyper-rectangle) if there's a part of further_hr within
	// sqrt( max_dist_sq ) of pt.
	if ( hr_intersects_hs( further_hr, pt, std::sqrt( max_dist_sq ) ) ) {
		current->distance( sq_vec_distance( current->location(), pt ) );
		if ( current->distance() < neighbors.worst_distance() ) {
			neighbors.insert( current->point() );
			max_dist_sq = neighbors.worst_distance();
		}
		//std::cout << "matching ";
		//print_point( std::cout, pt );
		//std::cout << " with ";
		//print_point( std::cout, current->location() );
		//std::cout << std::endl;
		//std::cout << "distance = " << current->distance() << std::endl;
		//std::cout << "worst_distance = " << neighbors.worst_distance() << std::endl;

		KDPointList temp_neighbors( neighbors.max_values() );
		temp_neighbors.distance_cutoff( neighbors.distance_cutoff() );

		//std::cout << "calling nearest_neighbors on further child" << std::endl;
		nearest_neighbors(
			further, pt, further_hr, max_dist_sq, temp_neighbors
		);

		neighbors.merge( temp_neighbors );
	}

} // nearest_neighbors

void
nearest_neighbor(
	KDNodeOP & current,
	utility::vector1< numeric::Real > const & pt,
	HyperRectangle & bounds,
	numeric::Real max_dist_sq,

	// returns:
	KDNodeOP & nearest,
	numeric::Real & dist_sq
) {
	using numeric::Size;
	using numeric::Real;
	using utility::vector1;

	if ( !current ) {
		dist_sq = REALLY_BIG_DISTANCE; // hopefully big enough for most spaces
		return;
	}
	Size split_axis  = current->split_axis();
	Real split_value = current->location()[split_axis];

	KDNodeOP nearer, further;
	// mid is a hyper-plane through this location and perpendicular
	// to the split axis
	vector1< Real >
		mid  ( current->location() ),
		lower( bounds.lower()      ),
		upper( bounds.upper()      );
	HyperRectangle left_hr ( lower, mid );
	HyperRectangle right_hr( upper, mid );

	HyperRectangle nearer_hr, further_hr;
	if ( pt[split_axis] <= split_value ) {
		nearer     = current->left_child();
		further    = current->right_child();
		nearer_hr  = left_hr;
		further_hr = right_hr;
	} else {
		nearer     = current->right_child();
		further    = current->left_child();
		nearer_hr  = left_hr;
		further_hr = right_hr;
	}

	nearest_neighbor( nearer, pt, nearer_hr, max_dist_sq, nearest, dist_sq );
	max_dist_sq = std::min( dist_sq, max_dist_sq );

	// we need to seach the right-hand child (furthest hyper-rectangle)
	// if there's a part of further_hr within sqrt( max_dist_sq ) of pt.
	if ( hr_intersects_hs( further_hr, pt, std::sqrt( max_dist_sq ) ) ) {
		numeric::Real this_dist_sq = sq_vec_distance( current->location(), pt );
		if ( this_dist_sq < dist_sq ) {
			nearest = current;
			dist_sq = this_dist_sq;
			max_dist_sq = dist_sq;
		}
		KDNodeOP temp_nearest;
		numeric::Real temp_dist_sq( REALLY_BIG_DISTANCE );
		nearest_neighbor(
			further, pt, further_hr, max_dist_sq, temp_nearest, temp_dist_sq
		);
		if ( temp_dist_sq < dist_sq ) {
			nearest = temp_nearest;
			dist_sq = temp_dist_sq;
		}
	} // hr_intersects_hs
} // nearest_neighbor

} // kdtree
} // numeric
