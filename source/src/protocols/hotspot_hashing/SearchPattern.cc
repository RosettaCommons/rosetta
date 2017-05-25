// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author
//
// Unit headers
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <protocols/hotspot_hashing/SearchPattern.hh>

namespace protocols {
namespace hotspot_hashing {

utility::vector1<core::kinematics::Stub> LSMSearchPattern::Searchpoints()
{
	utility::vector1<core::kinematics::Stub> searchpoints;

	Vector xunit = Vector(0, 0, 1);
	Matrix normal_rotation = rotation_matrix( lsmspec_.direction.cross(xunit), angle_of(lsmspec_.direction, xunit));

	for ( core::Real x = -max_radius_; x <= max_radius_; x += translocation_sampling_ ) {
		for ( core::Real y = -max_radius_; y <= max_radius_; y += translocation_sampling_ ) {
			if ( sqrt(x*x + y*y) <= max_radius_ ) {
				for ( core::Real z = 0; z <= max_distance_; z += distance_sampling_ ) {
					for ( core::Real angle = 0; angle < 360; angle += angle_sampling_ ) {
						Vector translation = Vector(x, y, z);
						Matrix rotation = numeric::x_rotation_matrix_degrees(angle);

						core::kinematics::Stub tp(rotation * normal_rotation, translation + lsmspec_.position);
						searchpoints.push_back(tp);
					}
				}
			}
		}
	}

	return searchpoints;
}

bool LSMSearchPattern::parse_lsm_spec(std::string lsmstring, VectorPair & lsmspec)
{
	std::vector< std::string > vectors;
	boost::split( vectors, lsmstring, boost::is_any_of(":") );

	if ( vectors.size() != 2 ) {
		return false;
	}

	std::vector< std::string > position;
	std::vector< std::string > direction;

	boost::split( position, vectors[0], boost::is_any_of(",") );
	boost::split( direction, vectors[1], boost::is_any_of(",") );

	if ( position.size() != 3 || direction.size() != 3 ) {
		return false;
	}

	using boost::lexical_cast;

	lsmspec = VectorPair(
		Vector(
		lexical_cast<core::Real>(position[0]),
		lexical_cast<core::Real>(position[1]),
		lexical_cast<core::Real>(position[2])),
		Vector(
		lexical_cast<core::Real>(direction[0]),
		lexical_cast<core::Real>(direction[1]),
		lexical_cast<core::Real>(direction[2])));

	return true;
}

utility::vector1<core::kinematics::Stub> RotationSearchPattern::generate_search_points()
{
	utility::vector1<core::kinematics::Stub> searchpoints;

	for ( core::Real x_sample = x_min; x_sample <= x_max; x_sample += x_sampling ) {
		for ( core::Real y_sample = y_min; y_sample <= y_max; y_sample += y_sampling ) {
			searchpoints.push_back(
				core::kinematics::Stub(
				numeric::x_rotation_matrix_degrees(x_sample) * numeric::y_rotation_matrix_degrees(y_sample),
				Vector(0)));
		}
	}

	return searchpoints;
}

utility::vector1<core::kinematics::Stub> SphericalRotationSearchPattern::generate_search_points()
{
	utility::vector1<core::kinematics::Stub> searchpoints;

	for ( core::Real polar_sample = polar_min; polar_sample <= polar_max; polar_sample += polar_rotation_sampling ) {
		for ( core::Real altitude_sample = altitude_min; altitude_sample <= altitude_max; altitude_sample += altitude_rotation_sampling ) {
			for ( core::Real azmiuth_sample = azmiuth_min; azmiuth_sample <= azmiuth_max; azmiuth_sample += azmiuth_rotation_sampling ) {
				searchpoints.push_back(
					core::kinematics::Stub(
					numeric::x_rotation_matrix_degrees(azmiuth_sample) * numeric::y_rotation_matrix_degrees(altitude_sample) * numeric::x_rotation_matrix_degrees(polar_sample),
					Vector(0)));
			}
		}
	}

	return searchpoints;
}

utility::vector1<core::kinematics::Stub> CartesianSearchPattern::generate_search_points()
{
	utility::vector1<core::kinematics::Stub> searchpoints;

	if ( x_sampling <= 0 ) {
		x_sampling = x_max - x_min + 1;
	}

	if ( y_sampling <= 0 ) {
		y_sampling = y_max - y_min + 1;
	}

	if ( z_sampling <= 0 ) {
		z_sampling = z_max - z_min + 1;
	}

	for ( core::Real x_sample = x_min; x_sample <= x_max; x_sample += x_sampling ) {
		for ( core::Real y_sample = y_min; y_sample <= y_max; y_sample += y_sampling ) {
			for ( core::Real z_sample = z_min; z_sample <= z_max; z_sample += z_sampling ) {
				searchpoints.push_back(
					core::kinematics::Stub(
					Matrix::identity(),
					Vector(x_sample, y_sample, z_sample)));
			}
		}
	}

	return searchpoints;
}

utility::vector1<core::kinematics::Stub> PartitionedSearchPattern::Searchpoints()
{
	utility::vector1<core::kinematics::Stub> sourcepoints = source_pattern_->Searchpoints();

	utility::vector1<core::kinematics::Stub> searchpoints;
	searchpoints.reserve((sourcepoints.size() / total_partitions_) + 1);

	for ( core::Size i = partition_; i < sourcepoints.size(); i += total_partitions_ ) {
		searchpoints.push_back(sourcepoints[i+1]);
	}

	return searchpoints;
}

utility::vector1<core::kinematics::Stub> ComposeSearchPatterns::Searchpoints()
{
	utility::vector1<core::kinematics::Stub> sourcepoints_a = source_pattern_a->Searchpoints();
	utility::vector1<core::kinematics::Stub> sourcepoints_b = source_pattern_b->Searchpoints();

	utility::vector1<core::kinematics::Stub> searchpoints;
	searchpoints.reserve(sourcepoints_a.size() * sourcepoints_b.size());

	for ( core::Size a = 1; a <= sourcepoints_a.size(); a += 1 ) {
		for ( core::Size b = 1; b <= sourcepoints_b.size(); b += 1 ) {
			core::kinematics::Stub result;
			core::kinematics::RT(core::kinematics::default_stub, sourcepoints_b[b]).make_jump(sourcepoints_a[a], result);
			searchpoints.push_back(result);
		}
	}

	return searchpoints;
}

utility::vector1<core::kinematics::Stub> SearchPatternTransform::Searchpoints()
{
	utility::vector1<core::kinematics::Stub> sourcepoints = source_pattern_->Searchpoints();

	utility::vector1<core::kinematics::Stub> searchpoints;
	searchpoints.reserve(sourcepoints.size());

	for ( core::Size i = 1; i <= sourcepoints.size(); i++ ) {
		searchpoints.push_back( transform(sourcepoints[i]));
	}

	return searchpoints;
}

utility::vector1<core::kinematics::Stub> SearchPatternExpansion::Searchpoints()
{
	utility::vector1<core::kinematics::Stub> sourcepoints = source_pattern_->Searchpoints();

	utility::vector1<core::kinematics::Stub> searchpoints;

	for ( core::Size i = 1; i <= sourcepoints.size(); i++ ) {
		utility::vector1<core::kinematics::Stub> newpoints = expand(sourcepoints[i]);

		searchpoints.insert(searchpoints.end(), newpoints.begin(), newpoints.end());
	}

	return searchpoints;
}

}
}

