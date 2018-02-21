// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/KMedoidsOnTheFly.hh
/// @brief Clustering algorithm used in multistage_rosetta_scripts.
/// @details This algorithm is stochastic but usually pretty good at
/// finding good clusters quickly. The methods in this file do not precalculate any distances.
/// Instead, all distances are measured "on the fly". In practice, this does not appear to
/// be a measurable slowdown.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_multistage_rosetta_scripts_cluster_KMedoidsOnTheFly_HH
#define INCLUDED_protocols_multistage_rosetta_scripts_cluster_KMedoidsOnTheFly_HH

#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetric.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {

///@brief Finds "num_medoids" number of clusters from within the "points" vector.
///If a point becomes the center of a cluster, it's index in "medoid_positions" is set to true.
///False otherwise. medoid_positions.size() will equal points.size()
void
k_medoids_on_the_fly (
	utility::vector1< ClusterMetricCOP > const & points,
	unsigned short int const num_medoids, //A.K.A. "k"
	utility::vector1< bool > & medoid_positions
);

///@brief This overload provides you with the vector of results.
/// It is mildly more convenient but mildly less efficient.
/// Differences are mild all around.
inline utility::vector1< bool >
k_medoids_on_the_fly (
	utility::vector1< ClusterMetricCOP > const & points,
	unsigned short int const num_medoids //A.K.A. "k"
) {
	utility::vector1< bool > medoid_positions;
	k_medoids_on_the_fly( points, num_medoids, medoid_positions );
	return medoid_positions;
}

} // cluster
} // multistage_rosetta_scripts
} // utility

#endif
