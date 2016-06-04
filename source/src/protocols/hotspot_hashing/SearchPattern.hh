// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_protocols_hotspot_hashing_SearchPattern_hh
#define INCLUDED_protocols_hotspot_hashing_SearchPattern_hh

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/exit.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#include <core/types.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>

#include <protocols/hotspot_hashing/SearchPattern.fwd.hh>

namespace protocols {
namespace hotspot_hashing {

typedef numeric::xyzMatrix< core::Real > Matrix;
typedef numeric::xyzVector< core::Real > Vector;


class VectorPair
{
public:
	VectorPair() :
		position(),
		direction()
	{
	};

	VectorPair(Vector position, Vector direction) :
		position(position),
		direction(direction)
	{
	};

	VectorPair(const VectorPair& src) :
		position(src.position),
		direction(src.direction)
	{
	};

	Vector position;
	Vector direction;
};

inline std::ostream& operator<<(std::ostream &stream, const Vector &vector)
{
	return stream <<  vector.x() << "," << vector.y() << "," << vector.z();
}

inline std::istream &
operator >>( std::istream & stream, Vector & v )
{
	core::Real x, y, z;
	stream >> x;
	stream.ignore(1, ',');
	stream >> y;
	stream.ignore(1, ',');
	stream >> z;

	v.x(x);
	v.y(y);
	v.z(z);

	return stream;
}

inline std::ostream& stub_to_points(std::ostream &stream, const core::kinematics::Stub &stub)
{
	/// Generate three points that would yield this stub
	Vector a = stub.local2global(Vector(0, 0, 0));
	Vector b = stub.local2global(Vector(-1, 0, 0));
	Vector c = stub.local2global(Vector(0, 1, 0));

	return stream << a << ";" << b << ";" << c;
}

inline std::istream& stub_from_points( std::istream & stream, core::kinematics::Stub & stub )
{
	Vector a, b, c;
	stream >> a;
	stream.ignore(1, ';');
	stream >> b;
	stream.ignore(1, ';');
	stream >> c;

	stub.from_four_points(a, a, b, c);

	return stream;
}

inline std::ostream& operator<<(std::ostream &stream, const VectorPair &pair)
{
	return stream << pair.position << ";" << pair.direction;
}

inline std::istream &
operator >>( std::istream & stream, VectorPair & pair )
{
	stream >> pair.position;
	stream.ignore(1, ';');
	stream >> pair.direction;

	return stream;
}

class SearchPattern : public utility::pointer::ReferenceCount
{
public:
	virtual utility::vector1<core::kinematics::Stub> Searchpoints() = 0;
};

class ConstPattern : public SearchPattern
{
public:
	ConstPattern() :
		target_stub(core::kinematics::default_stub)
	{}

	ConstPattern(core::kinematics::Stub target) :
		target_stub(target)
	{}

	virtual utility::vector1<core::kinematics::Stub> Searchpoints()
	{
		utility::vector1<core::kinematics::Stub> searchpoints;
		searchpoints.push_back(target_stub);
		return searchpoints;
	}

	core::kinematics::Stub target_stub;
};

class TestPattern : public SearchPattern
{
public:
	TestPattern()
	{
	}

	virtual utility::vector1<core::kinematics::Stub> Searchpoints()
	{
		utility::vector1<core::kinematics::Stub> searchpoints;

		core::Real x = 0;
		core::Real y = 0;
		core::Real z = 0;

		for ( x = -1; x <= 1; x += 1 ) {
			for ( y = -1; y <= 1; y += 1 ) {
				for ( z = -1; z <= 1; z += 1 ) {
					for ( core::Real angle = 0; angle < 360; angle += 180 ) {
						Vector translation = Vector(x, y, z);
						Matrix rotation = numeric::z_rotation_matrix_degrees(angle);

						core::kinematics::Stub tp(rotation, translation);
						searchpoints.push_back(tp);
					}
				}
			}
		}

		return searchpoints;
	};
};

class LSMSearchPattern : public SearchPattern
{
public:
	LSMSearchPattern(
		VectorPair lsmspec,
		core::Real angle_sampling,
		core::Real translocation_sampling,
		core::Real max_radius,
		core::Real distance_sampling,
		core::Real max_distance) :
		lsmspec_(lsmspec),
		angle_sampling_(angle_sampling),
		translocation_sampling_(translocation_sampling),
		max_radius_(max_radius),
		distance_sampling_(distance_sampling),
		max_distance_(max_distance)
	{

	}

	virtual utility::vector1<core::kinematics::Stub> Searchpoints()
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

	static bool parse_lsm_spec(std::string lsmstring, VectorPair & lsmspec)
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


private:
	VectorPair lsmspec_;
	core::Real angle_sampling_;
	core::Real translocation_sampling_;
	core::Real max_radius_;
	core::Real distance_sampling_;
	core::Real max_distance_;
};

class RotationSearchPattern : public SearchPattern
{
public:
	RotationSearchPattern(
		core::Real x_sampling,
		core::Real y_sampling,
		core::Real x_min = 0,
		core::Real x_max = 360,
		core::Real y_min = 0,
		core::Real y_max = 90) :
		x_sampling(x_sampling),
		y_sampling(y_sampling),
		x_min(x_min),
		x_max(x_max),
		y_min(y_min),
		y_max(y_max),
		searchpoints_()
	{
		searchpoints_ = generate_search_points();
	}

	core::Real x_sampling;
	core::Real y_sampling;
	core::Real x_min;
	core::Real x_max;
	core::Real y_min;
	core::Real y_max;

	utility::vector1<core::kinematics::Stub> generate_search_points()
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

	utility::vector1<core::kinematics::Stub> Searchpoints()
	{
		return searchpoints_;
	}

	utility::vector1<core::kinematics::Stub> searchpoints_;
};

class SphericalRotationSearchPattern : public SearchPattern
{
	// http://upload.wikimedia.org/wikipedia/commons/4/4f/3D_Spherical.svg
	// Uses spherical coordinates as outlined in the above diagram.
	// Altitude rotate is equiv. to theta, rotation off the polar axis.
	// Azmiuth rotate is equiv. to phi, rotation around the polar axis.
	// Polar rotate is rotation around the result vector from alt/azm rotation.

public:
	SphericalRotationSearchPattern(
		core::Real polar_rotation_sampling,
		core::Real altitude_rotation_sampling,
		core::Real azmiuth_rotation_sampling,
		core::Real polar_min = 0,
		core::Real polar_max = 360,
		core::Real altitude_min = 0,
		core::Real altitude_max = 360,
		core::Real azmiuth_min = 0,
		core::Real azmiuth_max = 360) :
		polar_rotation_sampling(polar_rotation_sampling),
		altitude_rotation_sampling(altitude_rotation_sampling),
		azmiuth_rotation_sampling(azmiuth_rotation_sampling),
		polar_min(polar_min),
		polar_max(polar_max),
		altitude_min(altitude_min),
		altitude_max(altitude_max),
		azmiuth_min(azmiuth_min),
		azmiuth_max(azmiuth_max),
		searchpoints_()
	{
		searchpoints_ = generate_search_points();
	}

	core::Real polar_rotation_sampling;
	core::Real altitude_rotation_sampling;
	core::Real azmiuth_rotation_sampling;
	core::Real polar_min;
	core::Real polar_max;
	core::Real altitude_min;
	core::Real altitude_max;
	core::Real azmiuth_min;
	core::Real azmiuth_max;

	utility::vector1<core::kinematics::Stub> generate_search_points()
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

	utility::vector1<core::kinematics::Stub> Searchpoints()
	{
		return searchpoints_;
	}

	utility::vector1<core::kinematics::Stub> searchpoints_;
};

class CartesianSearchPattern : public SearchPattern
{
public:
	CartesianSearchPattern(
		core::Real x_sampling,
		core::Real y_sampling,
		core::Real z_sampling,
		core::Real x_min,
		core::Real x_max,
		core::Real y_min,
		core::Real y_max,
		core::Real z_min,
		core::Real z_max) :
		x_sampling(x_sampling),
		y_sampling(y_sampling),
		z_sampling(z_sampling),
		x_min(x_min),
		x_max(x_max),
		y_min(y_min),
		y_max(y_max),
		z_min(z_min),
		z_max(z_max),
		searchpoints_()
	{
		searchpoints_ = generate_search_points();
	}

	core::Real x_sampling;
	core::Real y_sampling;
	core::Real z_sampling;
	core::Real x_min;
	core::Real x_max;
	core::Real y_min;
	core::Real y_max;
	core::Real z_min;
	core::Real z_max;

	utility::vector1<core::kinematics::Stub> generate_search_points()
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

	utility::vector1<core::kinematics::Stub> Searchpoints()
	{
		return searchpoints_;
	}

	utility::vector1<core::kinematics::Stub> searchpoints_;
};

class PartitionedSearchPattern : public SearchPattern
{
public:
	PartitionedSearchPattern(
		SearchPatternOP source_pattern,
		core::Size partition,
		core::Size total_partitions) :
		source_pattern_(source_pattern),
		partition_(partition),
		total_partitions_(total_partitions)
	{
		runtime_assert(partition < total_partitions);
	}

	virtual utility::vector1<core::kinematics::Stub> Searchpoints()
	{
		utility::vector1<core::kinematics::Stub> sourcepoints = source_pattern_->Searchpoints();

		utility::vector1<core::kinematics::Stub> searchpoints;
		searchpoints.reserve((sourcepoints.size() / total_partitions_) + 1);

		for ( core::Size i = partition_; i < sourcepoints.size(); i += total_partitions_ ) {
			searchpoints.push_back(sourcepoints[i+1]);
		}

		return searchpoints;
	}

private:
	SearchPatternOP source_pattern_;
	core::Size partition_;
	core::Size total_partitions_;
};

class ComposeSearchPatterns : public SearchPattern
{
public:
	ComposeSearchPatterns(
		SearchPatternOP source_pattern_a,
		SearchPatternOP source_pattern_b) :
		source_pattern_a(source_pattern_a),
		source_pattern_b(source_pattern_b)
	{ }

	virtual utility::vector1<core::kinematics::Stub> Searchpoints()
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


private:
	SearchPatternOP source_pattern_a;
	SearchPatternOP source_pattern_b;

};

class SearchPatternTransform : public SearchPattern
{
public:
	SearchPatternTransform( SearchPatternOP source_pattern ) :
		source_pattern_(source_pattern)
	{}

	virtual utility::vector1<core::kinematics::Stub> Searchpoints()
	{
		utility::vector1<core::kinematics::Stub> sourcepoints = source_pattern_->Searchpoints();

		utility::vector1<core::kinematics::Stub> searchpoints;
		searchpoints.reserve(sourcepoints.size());

		for ( core::Size i = 1; i <= sourcepoints.size(); i++ ) {
			searchpoints.push_back( transform(sourcepoints[i]));
		}

		return searchpoints;
	}

	virtual core::kinematics::Stub transform(core::kinematics::Stub source) = 0;

private:
	SearchPatternOP source_pattern_;
};

class SearchPatternExpansion : public SearchPattern
{
public:
	SearchPatternExpansion( SearchPatternOP source_pattern ) :
		source_pattern_(source_pattern)
	{}

	virtual utility::vector1<core::kinematics::Stub> Searchpoints()
	{
		utility::vector1<core::kinematics::Stub> sourcepoints = source_pattern_->Searchpoints();

		utility::vector1<core::kinematics::Stub> searchpoints;

		for ( core::Size i = 1; i <= sourcepoints.size(); i++ ) {
			utility::vector1<core::kinematics::Stub> newpoints = expand(sourcepoints[i]);

			searchpoints.insert(searchpoints.end(), newpoints.begin(), newpoints.end());
		}

		return searchpoints;
	}

	virtual utility::vector1<core::kinematics::Stub> expand(core::kinematics::Stub source) = 0;

private:
	SearchPatternOP source_pattern_;
};

}
}

#endif
