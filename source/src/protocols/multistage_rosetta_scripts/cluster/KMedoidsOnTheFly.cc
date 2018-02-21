// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/KMedoidsOnTheFly.cc
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <protocols/multistage_rosetta_scripts/cluster/KMedoidsOnTheFly.hh>
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

using uint = unsigned int;

///@brief This assumes you already created k "true"s in medoids vector.
/// Just goes through and updates medoid_for_cluster and cluster_for_point (but only for medoids)
void
assign_initial_medoids_on_the_fly(
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
void
assign_medoids_on_the_fly (
	utility::vector1< ClusterMetricCOP > const & points,
	utility::vector1< bool > & medoids,
	utility::vector1< uint > & medoid_for_cluster,
	utility::vector1< uint > & cluster_for_point
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
			cost += points[ point ]->distance( * points[ other_point ] );
		}
		if ( new_medoids[ cluster ].first == 0 || new_medoids[ cluster ].second > cost ) {
			new_medoids[ cluster ].first = point;
			new_medoids[ cluster ].second = cost;
		}
	}

	//reset:
	medoids.assign( medoids.size(), false );

	for ( uint cluster = 1; cluster <= num_clusters; ++cluster ) {
		uint const medoid = new_medoids[ cluster ].first;
		medoids[ medoid ] = true;
		cluster_for_point[ medoid ] = cluster;
		medoid_for_cluster[ cluster ] = medoid;
	}
}

///@brief After the medoids have been assigned, go through each point and determine which cluster it is in.
/// Returns total distance cost.
double
assign_non_medoids_to_clusters_on_the_fly (
	utility::vector1< ClusterMetricCOP > const & points,
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
			uint const other_point = medoid_for_cluster[ cluster ];
			double const dist = points[ i ]->distance( * points[ other_point ] );

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
k_medoids_on_the_fly (
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

	//Determine initial medoids:
	utility::vector1< uint > medoid_for_cluster( num_medoids, 0 );
	utility::vector1< uint > cluster_for_point( num_points, 0 );
	assign_initial_medoids_on_the_fly( medoids, medoid_for_cluster, cluster_for_point );
	double current_cost =
		assign_non_medoids_to_clusters_on_the_fly( points, medoids, medoid_for_cluster, cluster_for_point );

	best_medoids = medoids;
	debug_assert( vector_has_n_trues( best_medoids, num_medoids ) );
	debug_assert( vector_has_n_trues( medoids, num_medoids ) );

	while ( true ) {
		assign_medoids_on_the_fly( points, medoids, medoid_for_cluster, cluster_for_point );
		double const cost =
			assign_non_medoids_to_clusters_on_the_fly( points, medoids, medoid_for_cluster, cluster_for_point );

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
