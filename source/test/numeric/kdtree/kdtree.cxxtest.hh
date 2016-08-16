// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/geometry/kdtree.cxxtest.hh
/// @brief  test suite for kd-tree code
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <utility/vector1.hh>

#include <numeric/types.hh>

#include <numeric/kdtree/util.hh>
#include <numeric/kdtree/calc_distance.hh>
#include <numeric/kdtree/construct_kdtree.hh>
#include <numeric/kdtree/nearest_neighbors.hh>
#include <numeric/kdtree/KDTree.hh>
#include <numeric/kdtree/KDNode.hh>
#include <numeric/kdtree/KDPointList.hh>
#include <numeric/kdtree/KDNode.fwd.hh>
#include <numeric/kdtree/HyperRectangle.hh>
#include <numeric/kdtree/WrappedReal.hh>

// --------------- Test Class --------------- //

class KDTreeTests : public CxxTest::TestSuite {

public:

	utility::vector1< numeric::Real > query;
	utility::vector1< utility::vector1< numeric::Real > > points;
	utility::vector1< utility::pointer::ReferenceCountOP > data;

	utility::vector1< numeric::kdtree::KDPointOP > heap_points;

	// shared initialization
	void setUp() {
		using numeric::Real;
		using utility::vector1;
		using namespace numeric::kdtree;

		static bool init( false );
		if ( !init ) {
			query.resize( 2, 0.0 );
			points.resize( 8, query );
			points[1][1] = 1;
			points[1][2] = 7; // 50
			points[2][1] = 8;
			points[2][2] = 4; // 80
			points[3][1] = 3;
			points[3][2] = 3; // 18
			points[4][1] = 7;
			points[4][2] = 8; // 113
			points[5][1] = 3;
			points[5][2] = 4; // 25
			points[6][1] = 6;
			points[6][2] = 7; // 85
			points[7][1] = 1;
			points[7][2] = 8; // 65
			points[8][1] = 1;
			points[8][2] = 1; // 2

			data.resize( points.size(), NULL );
			for ( numeric::Size ii = 1; ii <= points.size(); ++ii ) {
				WrappedRealOP val( new WrappedReal( ii ) );
				data[ii] = val;
			}

			init = true;
		} // if ( !init )
	} // setUp

	// shared finalization
	void tearDown() {
	}

	// --------------- Test Cases --------------- //

	/// @brief test for naive nearest-neighbors search
	void test_naive_nn() {
		using numeric::Real;
		using utility::vector1;
		using namespace numeric::kdtree;

		typedef vector1< vector1< Real > >::const_iterator iter;
		vector1< Real > sq_distances;
		for ( iter it = points.begin(), end = points.end(); it != end; ++it ) {
			sq_distances.push_back( sq_vec_distance( *it, query ) );
		}

		TS_ASSERT( sq_distances[1] == 50 );
		TS_ASSERT( sq_distances[2] == 80 );
		TS_ASSERT( sq_distances[3] == 18 );
		TS_ASSERT( sq_distances[4] == 113 );
		TS_ASSERT( sq_distances[5] == 25 );
		TS_ASSERT( sq_distances[6] == 85 );
		TS_ASSERT( sq_distances[7] == 65 );
		TS_ASSERT( sq_distances[8] == 2 );
	}

	/// @brief test tree construction and simple nearest-neighbor
	/// search.
	void test_basic_kdtree() {
		using numeric::Real;
		using utility::vector1;
		using namespace numeric::kdtree;

		KDTree tree( points, data );
		KDNodeOP nearest;
		numeric::Real dist_sq;

		TS_ASSERT( tree.size() == 8 );
		TS_ASSERT( tree.ndim() == 2 );

		nearest_neighbor( tree, query, nearest, dist_sq );
		vector1< Real > closest( 2, 1 );

		TS_ASSERT( dist_sq == 2 );
		TS_ASSERT( nearest->location() == closest );
		WrappedRealOP val = utility::pointer::dynamic_pointer_cast< WrappedReal > ( nearest->data() );
		TS_ASSERT( val->val() == 8.0 );
	}

	/// @brief test KDPointList operations.
	void test_points() {
		using numeric::Real;
		using numeric::Size;
		using utility::vector1;
		using namespace numeric::kdtree;

		Size const ndim( 4 );
		Size const wanted( 5 );
		KDPointList container( wanted );
		vector1< Real > origin( ndim, 0.0 );

		// add some randomly high values for consideration.
		container.insert(
			KDPointOP( new KDPoint( vector1< Real >( 4, 25 ),
			sq_vec_distance( vector1< Real >( 4, 25 ), origin ) ) )
		);
		container.insert(
			KDPointOP( new KDPoint( vector1< Real >( 4, 67 ),
			sq_vec_distance( vector1< Real >( 4, 67 ), origin ) ) )
		);
		container.insert(
			KDPointOP( new KDPoint( vector1< Real >( 4, 99 ),
			sq_vec_distance( vector1< Real >( 4, 99 ), origin ) ) )
		);
		for ( Size ii = 1; ii <= 20; ++ii ) {
			vector1< Real > k_loc( ndim, ii );
			KDPointOP pt( new KDPoint( k_loc, sq_vec_distance( k_loc, origin ) ) );

			container.insert( pt );
		}

		TS_ASSERT( container.size() == wanted );
		TS_ASSERT( container.worst_distance() == 100 );

		utility::vector1< KDPointOP > values = container.sorted_values();
		//TS_ASSERT( container.worst().location() == vector1< Real >( ndim, 5 ) );
		for ( Size ii = 1; ii <= wanted; ++ii ) {
			vector1< Real > k_loc( ndim, ii );
			KDPoint pt( k_loc, sq_vec_distance( k_loc, origin ) );

			TS_ASSERT( *values[ii] == pt );
		}

		// make container2 with one value that would fit into container,
		// and three that would not.
		{ // scope for merging test
			KDPointList container2( wanted );
			container2.insert(
				KDPointOP( new KDPoint( vector1< Real >( 4, 2.5 ),
				sq_vec_distance( vector1< Real >( 4, 2.5 ), origin ) ) )
			);
			container2.insert(
				KDPointOP( new KDPoint( vector1< Real >( 4, 25 ),
				sq_vec_distance( vector1< Real >( 4, 25 ), origin ) ) )
			);
			container2.insert(
				KDPointOP( new KDPoint( vector1< Real >( 4, 67 ),
				sq_vec_distance( vector1< Real >( 4, 67 ), origin ) ) )
			);
			container2.insert(
				KDPointOP( new KDPoint( vector1< Real >( 4, 99 ),
				sq_vec_distance( vector1< Real >( 4, 99 ), origin ) ) )
			);
			container2.insert(
				KDPointOP( new KDPoint( vector1< Real >( 4, 99 ),
				sq_vec_distance( vector1< Real >( 4, 99 ), origin ) ) )
			);
			TS_ASSERT( container2.size() == 5 );
			TS_ASSERT( container2.worst_distance() == 39204 );
			// test merging
			container.merge( container2 );


			TS_ASSERT( container.worst_distance() == 64 );

			values = container.sorted_values();
			TS_ASSERT( values.size() == wanted );
			TS_ASSERT( values[1]->distance() ==  4 );
			TS_ASSERT( values[2]->distance() == 16 );
			TS_ASSERT( values[3]->distance() == 25 );
			TS_ASSERT( values[4]->distance() == 36 );
			TS_ASSERT( values[5]->distance() == 64 );
		} // scope for merging test

	} // test_points

	/// test case for nearest_neighbors search using KDPointList
	void test_nearest_neighbors() {
		using numeric::Real;
		using numeric::Size;
		using utility::vector1;
		using namespace numeric::kdtree;

		KDTree tree( points, data );

		TS_ASSERT( tree.size() == 8 );
		TS_ASSERT( tree.ndim() == 2 );

		Size const wanted( 5 );
		KDPointList neighbors = nearest_neighbors( tree, query, wanted );
		TS_ASSERT( neighbors.size() == wanted );

		vector1< KDPointOP > values = neighbors.sorted_values();
		TS_ASSERT( values[1]->location() == vector1< Real >( 2,1 ) );

		TS_ASSERT( values[1]->distance() == 2  );
		TS_ASSERT( values[2]->distance() == 18 );
		TS_ASSERT( values[3]->distance() == 25 );
		TS_ASSERT( values[4]->distance() == 50 );
		TS_ASSERT( values[5]->distance() == 65 );

		WrappedRealOP nearest_val = utility::pointer::dynamic_pointer_cast< WrappedReal > ( values[1]->data() );
		TS_ASSERT( nearest_val->val() == 8.0 );
	}

	// test upper distance cutoff in NN search
	void test_nn_dist() {
		using numeric::Real;
		using numeric::Size;
		using utility::vector1;
		using namespace numeric::kdtree;

		KDTree tree( points, data );

		TS_ASSERT( tree.size() == 8 );
		TS_ASSERT( tree.ndim() == 2 );

		Size const wanted( 5 );
		KDPointList neighbors = nearest_neighbors( tree, query, wanted );
		TS_ASSERT( neighbors.size() == wanted );

		vector1< KDPointOP > values = neighbors.sorted_values();
		TS_ASSERT( values[1]->location() == vector1< Real >( 2,1 ) );

		TS_ASSERT( values[1]->distance() == 2  );
		TS_ASSERT( values[2]->distance() == 18 );
		TS_ASSERT( values[3]->distance() == 25 );
		TS_ASSERT( values[4]->distance() == 50 );
		TS_ASSERT( values[5]->distance() == 65 );

		WrappedRealOP nearest_val = utility::pointer::dynamic_pointer_cast< WrappedReal > ( values[1]->data() );
		TS_ASSERT( nearest_val->val() == 8.0 );
	}

	// linear transformation that maps values into range (0,1)
	void test_linear_transformation() {
		using numeric::Real;
		using numeric::Size;
		using utility::vector1;
		using namespace numeric::kdtree;

		vector1< vector1< Real > > values;

		for ( Size ii = 10; ii <= 90; ++ii ) {
			values.push_back( vector1< Real >( ii, ii + 1 ) );
		}

		transform_percentile( values );
		TS_ASSERT( values.front()[1] == 0 );
		TS_ASSERT( values.front()[2] == 0 );
		TS_ASSERT( values.back() [1] == 1 );
		TS_ASSERT( values.back() [2] == 1 );
		// maybe do some finer-grained testing here.
	}
}; // KDTreeTests
