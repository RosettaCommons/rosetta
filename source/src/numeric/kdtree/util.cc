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

#include <numeric/kdtree/util.hh>
#include <numeric/kdtree/KDNode.hh>
#include <numeric/kdtree/KDPoint.hh>
#include <numeric/kdtree/HyperRectangle.hh>
#include <numeric/kdtree/HyperRectangle.fwd.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>

#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

namespace numeric {
namespace kdtree {

HyperRectangleOP get_percentile_bounds(
	utility::vector1< utility::vector1< numeric::Real > > & points
) {
	using std::min;
	using std::max;
	using numeric::Real;
	using utility::vector1;

	// define lower/upper values
	vector1< Real >
		lower( points.front().begin(), points.front().end() ),
		upper( points.front().begin(), points.front().end() );

	for ( auto & point : points ) {
		for ( auto p_it = point.begin(), p_end = point.end(),
				l_it = lower.begin(), l_end = lower.end(),
				u_it = upper.begin(), u_end = upper.end();
				p_it != p_end && l_it != l_end && u_it != u_end;
				++p_it, ++l_it, ++u_it
				) {
			*l_it = std::min( *l_it, *p_it );
			*u_it = std::max( *u_it, *p_it );
		}
	} // rows

	return HyperRectangleOP( new HyperRectangle(
		upper, lower
		) );
}

void transform_percentile_single_pt(
	utility::vector1< numeric::Real > & point,
	HyperRectangleOP bounds
) {
	using numeric::Real;
	using utility::vector1;
	vector1< Real >
		lower( bounds->lower() ),
		upper( bounds->upper() );
	for ( auto p_it = point.begin(), p_end = point.end(),
			l_it = lower.begin(), l_end = lower.end(),
			u_it = upper.begin(), u_end = upper.end();
			p_it != p_end && l_it != l_end && u_it != u_end;
			++p_it, ++l_it, ++u_it
			) {
		*p_it = ( *p_it - *l_it ) / ( *u_it - *l_it );
	}
}

void transform_percentile(
	utility::vector1< utility::vector1< numeric::Real > > & points,
	HyperRectangleOP bounds
) {
	using std::min;
	using std::max;
	using numeric::Real;
	using utility::vector1;

	// define lower/upper values
	vector1< Real >
		lower( bounds->lower() ),
		upper( bounds->upper() );

	// transform values
	for ( auto & point : points ) {
		for ( auto p_it = point.begin(), p_end = point.end(),
				l_it = lower.begin(), l_end = lower.end(),
				u_it = upper.begin(), u_end = upper.end();
				p_it != p_end && l_it != l_end && u_it != u_end;
				++p_it, ++l_it, ++u_it
				) {
			*p_it = ( *p_it - *l_it ) / ( *u_it - *l_it );
		}
	} // rows
} // transform_percentile


void transform_percentile(
	utility::vector1< utility::vector1< numeric::Real > > & points
) {
	HyperRectangleOP bounds = get_percentile_bounds( points );
	transform_percentile( points, bounds );
}

utility::vector1< KDPointOP > make_points(
	utility::vector1< utility::vector1< numeric::Real > > const & points
) {

	using numeric::Real;
	using utility::vector1;

	vector1< KDPointOP > new_data;
	 for ( auto const & point : points ) {
		KDPointOP pt( new KDPoint( point ) );
		new_data.push_back( pt );
	}

	return new_data;
} // make_points

utility::vector1< KDPointOP > make_points(
	utility::vector1< utility::vector1< numeric::Real > > const & points,
	utility::vector1< utility::pointer::ReferenceCountOP > const & data
) {

	assert( points.size() == data.size() );

	using numeric::Real;
	using utility::vector1;
	using utility::pointer::ReferenceCountOP;

	vector1< KDPointOP > new_data;
	auto d_it = data.begin(), d_end = data.end();
	for ( auto p_it = points.begin(), p_end = points.end();
			p_it != p_end && d_it != d_end;
			++p_it, ++d_it
			) {
		KDPointOP pt( new KDPoint( *p_it, *d_it ) );
		new_data.push_back( pt );
	}

	return new_data;
} // make_points

void print_points(
	std::ostream & out,
	utility::vector1< utility::vector1< numeric::Real > > const & points
) {
	using numeric::Real;
	using utility::vector1;
	 for ( auto const & point : points ) {
		//for ( vector1< Real >::const_iterator val = pt->begin(),
		//  val_end = pt->end(); val != val_end; ++val
		//) {
		// out << ' ' << *val;
		//}
		print_point( out, point );
		out << std::endl;
	} // for points
}

void print_point(
	std::ostream & out,
	utility::vector1< numeric::Real > const & point
) {
	using numeric::Real;
	using utility::vector1;
	for ( double val : point ) {
		out << ' ' << val;
	}
}

// returns true if the given hyper-rectangle intersects with the given
// hypersphere.
bool hr_intersects_hs(
	HyperRectangle hr,
	utility::vector1< numeric::Real > const & pt,
	numeric::Real const r
) {
	using numeric::Size;
	using numeric::Real;
	using utility::vector1;

	vector1< Real >
		//qt( pt ),
		upper( hr.upper() ),
		lower( hr.lower() );

	Real dist_sq(0.0);
	//vector1< Real >::iterator qt_it = qt.begin();
	//std::cout << "testing intersection of sphere (";
	//print_point( std::cout, pt );
	//std::cout << "), r = " << r << std::endl;
	//std::cout << "lower = ";
	//print_point( std::cout, lower );
	//std::cout << std::endl;
	//std::cout << "upper = ";
	//print_point( std::cout, upper );
	//std::cout << std::endl;
	for ( vector1< Real >::const_iterator
			pt_it = pt.begin(), pt_end = pt.end(),
			lower_it = lower.begin(), lower_end = lower.end(),
			upper_it = upper.begin(), upper_end = upper.end();
			pt_it != pt_end && lower_it != lower_end && upper_it != upper_end;
			++pt_it, ++lower_it, ++upper_it //, qt_it
			) {
		using std::pow;

		//if ( *pt_it <= *lower_it || *pt_it >= *upper_it ) {
		// std::cout << "points are "
		//  << *pt_it << "," << *lower_it << "," << *upper_it
		//  << std::endl;
		// dist_sq += std::min(
		//  pow( *pt_it - *lower_it, 2 ), pow( *pt_it - *upper_it, 2 )
		// );
		// std::cout << "(choosing between "
		//  << pow( *pt_it - *lower_it, 2 ) << " and "
		//  << pow( *pt_it - *upper_it, 2 ) << ")"
		//  << std::endl;
		//}
		//std::cout << "dist_sq = " << dist_sq << std::endl;
	}

	//std::cout << "dist_sq = " << dist_sq << std::endl;
	//std::cout << "r^2 = " << r * r << std::endl;

	return ( dist_sq <= (r * r) );
} // hr_intersects_hs

void print_tree(
	std::ostream & out,
	KDNodeOP const & current,
	Size current_depth,
	Size const width
) {
	for ( Size ii = 0; ii < current_depth * width; ++ii ) {
		out << " ";
	}

	if ( !current ) {
		out << "empty" << std::endl;
		return;
	}

	out << "point ";
	print_point( out, current->location() );
	std::cout << " split_pt = " << current->split_axis() << ", "
		<< current->location()[ current->split_axis() ];
	out << std::endl;

	for ( Size ii = 0; ii < current_depth * width; ++ii ) {
		out << " ";
	}
	out << "left child" << std::endl;
	print_tree( out, current->left_child(), current_depth + 1, width );

	for ( Size ii = 0; ii < current_depth * width; ++ii ) {
		out << " ";
	}
	out << "right child" << std::endl;
	print_tree( out, current->right_child(), current_depth + 1, width );
} // print_tree

bool is_legal_kdtree_below_node(
	KDNodeOP const & current
) {
	if ( !current ) { return true; }

	Size const axis( current->split_axis() );
	Real const value( current->location()[axis] );

	return(
		is_legal_less_than( current->left_child(), axis, value ) &&
		is_legal_greater_than( current->right_child(), axis, value )
	);
}

bool is_legal_less_than(
	KDNodeOP const & current,
	Size const axis,
	Real const value
) {
	if ( !current ) { return true; }

	Real const this_val( current->location()[axis] );
	bool const this_node_good( current->location()[axis] <= value );
	if ( !this_node_good ) {
		std::cout << "Error with value: " << this_val  << " <= " << value
			<< std::endl;
	}

	bool children_are_good(
		is_legal_less_than( current->left_child(), axis, value ) &&
		is_legal_less_than( current->right_child(), axis, value )
	);

	return ( this_node_good && children_are_good );
}

bool is_legal_greater_than(
	KDNodeOP const & current,
	Size const axis,
	Real const value
) {
	if ( !current ) { return true; }

	Real const this_val( current->location()[axis] );
	bool const this_node_good( this_val > value );

	if ( !this_node_good ) {
		std::cout << "Error with value: " << this_val  << " > " << value
			<< std::endl;
	}

	bool children_are_good(
		is_legal_greater_than( current->left_child(), axis, value ) &&
		is_legal_greater_than( current->right_child(), axis, value )
	);
	return ( this_node_good && children_are_good );
}

} // kdtree
} // numeric
