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

#include <devel/init.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <core/types.hh>

namespace protocols {
namespace hotspot_hashing {

typedef numeric::xyzMatrix< core::Real > Matrix;
typedef numeric::xyzVector< core::Real > Vector;

class TransformPair
{
public:
	TransformPair (Vector translation, Matrix rotation) :
		translation(translation),
		rotation(rotation)
	{
	};

	TransformPair() :
		translation(0, 0, 0),
		rotation(Matrix::diag(1, 1, 1))
	{
	};

	Vector translation;
	Matrix rotation;
};

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

inline std::ostream& operator<<(std::ostream &strm, const Vector &vector) {
  return strm << "[" << vector.x() << "," << vector.y() << "," << vector.z() << "]";
}

inline std::ostream& operator<<(std::ostream &strm, const Matrix &matrix) {
  return strm << "[" << matrix.row_x() << "," << matrix.row_y() << "," << matrix.row_z() << "]";
}

inline std::ostream& operator<<(std::ostream &strm, const TransformPair &tp) {
  return strm << "(" << tp.translation << "," << tp.rotation << ")";
}

inline std::ostream& operator<<(std::ostream &strm, const VectorPair &vp) {
  return strm << "(" << vp.position << "," << vp.direction << ")";
}



class SearchPattern : public utility::pointer::ReferenceCount
{
	public:
		virtual utility::vector1<TransformPair> Searchpoints() = 0;
};

typedef utility::pointer::owning_ptr<SearchPattern>  SearchPatternOP;
typedef utility::pointer::owning_ptr<SearchPattern const>  SearchPatternCOP;

class ConstPattern : public SearchPattern
{
	public:
		virtual utility::vector1<TransformPair> Searchpoints()
		{
			utility::vector1<TransformPair> searchpoints;
			searchpoints.push_back(TransformPair(Vector(10, 0, 0), numeric::z_rotation_matrix_degrees((core::Real)90)));
			searchpoints.push_back(TransformPair(Vector(0, 0, 0), numeric::z_rotation_matrix_degrees((core::Real)0)));
			searchpoints.push_back(TransformPair(Vector(-10, 0, 0), numeric::z_rotation_matrix_degrees((core::Real)270)));

			return searchpoints;
		}
};

class TestPattern : public SearchPattern
{
	public:
		TestPattern()
		{
		}

		virtual utility::vector1<TransformPair> Searchpoints()
		{
			utility::vector1<TransformPair> searchpoints;

			core::Real x = 0;
			core::Real y = 0;
			core::Real z = 0;

			for (x = -1; x <= 1; x += 1)
			{
				for (y = -1; y <= 1; y += 1)
				{
					for (z = -1; z <= 1; z += 1 )
					{
						for (core::Real angle = 0; angle < 360; angle += 180)
						{
							Vector translation = Vector(x, y, z);
							Matrix rotation = numeric::z_rotation_matrix_degrees(angle);

							TransformPair tp(translation, rotation);
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

		virtual utility::vector1<TransformPair> Searchpoints()
		{
			utility::vector1<TransformPair> searchpoints;

			Vector xunit = Vector(0, 0, 1);
			Matrix normal_rotation = rotation_matrix( lsmspec_.direction.cross(xunit), angle_of(lsmspec_.direction, xunit));
			
			for (core::Real x = -max_radius_; x <= max_radius_; x += translocation_sampling_)
			{
				for (core::Real y = -max_radius_; y <= max_radius_; y += translocation_sampling_)
				{
					if(sqrt(x*x + y*y) <= max_radius_)
					{
						for (core::Real z = 0; z <= max_distance_; z += distance_sampling_)
						{
							for (core::Real angle = 0; angle < 360; angle += angle_sampling_)
							{
								Vector translation = Vector(x, y, z);
								Matrix rotation = numeric::x_rotation_matrix_degrees(angle);

								TransformPair tp(translation + lsmspec_.position, rotation * normal_rotation);
								searchpoints.push_back(tp);
							}
						}
					}
				}
			}

			return searchpoints;
		}

	private:
		VectorPair lsmspec_;
		core::Real angle_sampling_;
		core::Real translocation_sampling_;
		core::Real max_radius_;
		core::Real distance_sampling_;
		core::Real max_distance_;
};


}
}

#endif
