// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/KMedoids.cc
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <protocols/multistage_rosetta_scripts/cluster/KMedoids.hh>
#include <algorithm>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
//#include <utility/random/random.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/testing.OptionKeys.gen.hh>

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {


/*
These methods depend on the precalculation of distances.
The distances we need to store can be organized in a triangle matrix
like the one shown in (1) below.

2D vectors can be tricky and wasteful, especially when they are not rectangular.
So I decided to get cute and store the info linearly as shown in (2) below.

To know how to access each element, I keep an array full of offsets for each row (3).
So now to access the element in row i, column j we look up triangle[ offset[ i ] + ( j - i ) ].

Example: If we want to look up element F (row 2, column 4), we get
triangle[ offset[ i ] + ( j - i ) ]
triangle[ offset[ 2 ] + ( 4 - 2 ) ]
triangle[ offset[ 2 ] + 2 ]
triangle[ 4 + 2 ]
triangle[ 6 ] (which is correct as shown in (2) )



(1) Triangle Matrix of Distances Between 5 Points:

.....1....2....3....4....5

..1       A    B    C    D

..2            E    F    G

..3                 H    I

..4                      J

..5


(2) Linear Representation of Triangle Matrix:

Value: A  B  C  D  E  F  G  H  I  J
Index: 1  2  3  4  5  6  7  8  9  10
.      ^           ^        ^     ^
.      |           |        |     |
.      |           |        |     Beginning for row 4
.      |           |        Beginning for row 3
.      |           Beginning for row 2
.      Beginning for row 1

(3) Vector of Offsets:

Offset: 0  4  7  9
Index:  1  2  3  4
*/

using uint = unsigned int;

///@brief Looks up the distance between two points in the tricky triangle matrix I hacked together
float distance(
	utility::vector1< float > const & triangle_matrix,
	utility::vector1< uint > const & offsets,
	uint pointA,
	uint pointB
) {

	if ( pointA == pointB ) return 0;

	if ( pointA > pointB ) {
		auto const difference = pointA - pointB;
		auto const index = offsets[ pointB ] + difference;
		return triangle_matrix[ index ];
	} else {
		auto const difference = pointB - pointA;
		auto const index = offsets[ pointA ] + difference;
		return triangle_matrix[ index ];
	}
}

///@brief maybe this is lazy programming on my part.
/// I just want to count how many possible interacitons there are without using math.
/// Feel free to replace this with a simple equation.
unsigned long long int
number_of_elements_in_exclusive_upper_triangle(
	unsigned short int const num_points
) {
	//maybe TODO find a more mathematical way to calculate this
	unsigned long long int count = 0;
	for ( uint i = 0; i < num_points; ++i ) {
		count += i;
	}
	return count;
}

///@brief Initialize the offset vector identified in schematic (2) above
void
initialize_offsets(
	utility::vector1< uint > & offsets,
	uint const num_points
){
	offsets[ 1 ] = 0;
	uint width = num_points - 1;
	for ( unsigned short int point = 2; point < num_points; ++point ) {
		offsets[ point ] = offsets[ point - 1 ] + width;
		--width;
	}
}

///@brief This assumes you already created k "true"s in medoids vector.
/// Just goes through and updates medoid_for_cluster and cluster_for_point (but only for medoids)
void
assign_initial_medoids(
	utility::vector1< bool > const & medoids,
	utility::vector1< uint > & medoid_for_cluster,
	utility::vector1< uint > & cluster_for_point
) {
	uint index = 0;
	for ( uint i = 1; i <= medoids.size(); ++i ) {
		if ( medoids[ i ] ) {
			medoid_for_cluster[ ++index ] = i;
			cluster_for_point[ i ] = index;
		}
	}
}

///@brief Iterates though every point and tests to see if it is the best medoid for its cluster.
void assign_medoids (
	utility::vector1< float > const & triangle_matrix,
	utility::vector1< uint > const & offsets,
	utility::vector1< bool > & medoids,
	utility::vector1< uint > & medoid_for_cluster,
	utility::vector1< uint > const & cluster_for_point
) {
	uint const num_points = medoids.size();
	uint const num_clusters = medoid_for_cluster.size();

	using point_candidate = std::pair< uint, float >;//id, score
	utility::vector1< point_candidate > new_medoids( num_clusters, point_candidate( 0, 0.0 ) );

	utility::vector1< utility::vector1< uint > > points_in_cluster;
	points_in_cluster.resize( num_clusters );

	for ( uint point = 1; point <= num_points; ++point ) {
		points_in_cluster[ cluster_for_point[ point ] ].push_back( point );
	}

	//Give every point the opportunity to be the medoid for its cluster
	for ( uint point = 1; point <= num_points; ++point ) {
		uint const cluster = cluster_for_point[ point ];
		float cost = 0.0;
		for ( uint other_point : points_in_cluster[ cluster ] ) {
			cost += distance( triangle_matrix, offsets, point, other_point );
		}
		if ( new_medoids[ cluster ].first == 0 || new_medoids[ cluster ].second > cost ) {
			new_medoids[ cluster ].first = point;
			new_medoids[ cluster ].second = cost;
		}
	}

	//re-assign medoids
	medoids.assign( medoids.size(), false );
	for ( uint cluster = 1; cluster <= num_clusters; ++cluster ) {
		uint const medoid = new_medoids[ cluster ].first;
		medoids[ medoid ] = true;
		runtime_assert( cluster_for_point[ medoid ] == cluster );
		medoid_for_cluster[ cluster ] = medoid;
	}
}

///@brief After the medoids have been assigned, go through each point and determine which cluster it is in.
/// Returns total distance cost.
double
assign_non_medoids_to_clusters (
	utility::vector1< float > const & triangle_matrix,
	utility::vector1< uint > const & offsets,
	utility::vector1< bool > const & medoids,
	utility::vector1< uint > const & medoid_for_cluster,
	utility::vector1< uint > & cluster_for_point
){
	double score = 0;

	uint const num_points = medoids.size();
	uint const num_clusters = medoid_for_cluster.size();

	for ( uint i = 1; i <= num_points; ++i ) {
		if ( medoids[ i ] ) continue;

		uint closest_cluster = 0;
		double smallest_score = 999999.9;

		for ( uint cluster = 1; cluster <= num_clusters; ++cluster ) {
			double const dist =
				distance( triangle_matrix, offsets, i, medoid_for_cluster[ cluster ] );
			if ( dist < smallest_score || closest_cluster == 0 ) {
				closest_cluster = cluster;
				smallest_score = dist;
			}
		}

		runtime_assert( closest_cluster );

		cluster_for_point[ i ] = closest_cluster;
		score += smallest_score;
	}

	return score;
}

///@brief Finds "num_medoids" number of clusters from within the "points" vector.
///If a point becomes the center of a cluster, it's index in "medoid_positions" is set to true.
///False otherwise. medoid_positions.size() will equal points.size()
/// @details This algorithm is stochastic but usually pretty good at
/// finding good clusters quickly. The methods in this file precalculate
/// all of the possible distances between the points in an attempt to save time.
/// In practice, this does not save any time at all (tested with SequenceMetrics of length 250)
/// so I recommend looking at protocols/multistage_rosetta_scripts/cluster/KMedoidsOnTheFly.hh
void
k_medoids_with_edge_precalculation (
	utility::vector1< ClusterMetricCOP > const & points,
	unsigned short int const num_medoids, //A.K.A. "k"
	utility::vector1< bool > & best_medoids
) {

	auto const num_points = points.size();

	if ( num_medoids >= num_points ) {
		best_medoids = utility::vector1< bool > ( num_points, true );
	}

	utility::vector1< bool > medoids ( num_points, false );
	for ( int i=1; i <= num_medoids; ++i ) {
		medoids[ i ] = true;
	}
	debug_assert( vector_has_n_trues( medoids, num_medoids ) );

	if ( basic::options::option[ basic::options::OptionKeys::testing::INTEGRATION_TEST ]() ) {
		std::shuffle( medoids.begin(), medoids.end(), std::default_random_engine( 0 ) );
	} else {
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::shuffle( medoids.begin(), medoids.end(), std::default_random_engine(seed) );
	}

	debug_assert( vector_has_n_trues( medoids, num_medoids ) );

	auto const num_comparisons = number_of_elements_in_exclusive_upper_triangle( num_points );
	utility::vector1< float > triangle_matrix( num_comparisons, 0.0 );
	utility::vector1< uint > offsets( num_points - 1 );
	initialize_offsets( offsets, num_points );

	//Populate matrix:
	for ( uint i = 1; i < num_points; ++i ) {
		for ( uint j = i + 1; j <= num_points; ++j ) {
			auto const difference = j - i;
			auto const index = offsets[ i ] + difference;
			triangle_matrix[ index ] = std::abs( points[ i ]->distance( * points[ j ] ) );
		}
	}

	//Determine initial medoids:
	utility::vector1< uint > medoid_for_cluster( num_medoids, 0 );
	utility::vector1< uint > cluster_for_point( num_points, 0 );
	assign_initial_medoids( medoids, medoid_for_cluster, cluster_for_point );
	double current_cost =
		assign_non_medoids_to_clusters( triangle_matrix, offsets, medoids, medoid_for_cluster, cluster_for_point );

	best_medoids = medoids;
	debug_assert( vector_has_n_trues( best_medoids, num_medoids ) );
	debug_assert( vector_has_n_trues( medoids, num_medoids ) );

	while ( true ) {
		assign_medoids( triangle_matrix, offsets, medoids, medoid_for_cluster, cluster_for_point );
		double const cost =
			assign_non_medoids_to_clusters( triangle_matrix, offsets, medoids, medoid_for_cluster, cluster_for_point );

		if ( cost < current_cost ) {
			current_cost = cost;
			debug_assert( vector_has_n_trues( medoids, num_medoids ) );
			best_medoids = medoids;
		} else {
			break;
		}
	}

	debug_assert( vector_has_n_trues( best_medoids, num_medoids ) );
}

} // cluster
} // multistage_rosetta_scripts
} // utility
