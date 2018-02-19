// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/Stub.hh
/// @brief  Stub class
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_Stub_hh
#define INCLUDED_core_kinematics_Stub_hh


// Unit headers
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/RT.fwd.hh>

// Numeric headers
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

// C++ headers
#include <cmath>

#include <core/types.hh>

namespace core {
namespace kinematics {


/////////////////////////////////////////////////////////////////////////////
/// @brief Stub class -- an object of orthogonal coordinate frame
///
/// @details an orthogonal coord frame M (matrix) centered at point V (vector),
/// defined by four points, one is for the center and the other three for
/// calculating the orthogonal frame. For example, a stub can be derived from
/// a backbone triplet N-CA-C centered at CA.
///
/// See @ref atomtree_overview "AtomTree overview and concepts" for details.
class Stub
{

public: // Types

	typedef  numeric::xyzMatrix< Real > Matrix;
	typedef  numeric::xyzVector< Real > Vector;

public: // Creation

	/// @brief constructor -- sets to "default" stub
	inline
	Stub():
		M( Matrix::identity() ),
		v( 0.0 )
	{}

	/// @brief copy constructor
	inline
	Stub(
		Matrix const & M_in,
		Vector const & v_in
	):
		M( M_in ),
		v( v_in )
	{}

	/// constructor from RT object
	Stub( RT const & rt );


	/// @brief constructor by four points
	///
	/// @details first point is the center (V) and the rest three are used to
	/// construct the coord frame (M). see member function from_four_points(...)
	/// construct a stub centered at v_in, as would come from building
	/// c then b then a
	inline
	Stub(
		Vector const & center,
		Vector const & a,
		Vector const & b,
		Vector const & c
	)
	{
		this->from_four_points( center, a, b, c );
	}

	/// @brief constructor by three points
	///
	/// @details first point is the center (V) and all the three are used to
	/// construct the coord frame (M). see member functon from_four_points(...)
	/// construct a stub as would come from building
	/// c then b then a
	inline
	Stub(
		Vector const & a,
		Vector const & b,
		Vector const & c
	)
	{
		this->from_four_points( a, a, b, c );
	}

public: // Methods

	/// @brief build a stub from a center and other three points a, b, c
	void
	from_four_points(
		Vector const & center,
		Vector const & a,
		Vector const & b,
		Vector const & c
	);

	/// @brief check if the stub is orthogonal under the tolerance cutoff
	bool
	is_orthogonal( double const & tolerance ) const
	{
		Matrix delta( M * M.transposed() - Matrix::identity() );
		return ( ( delta.col_x().length() +
			delta.col_y().length() +
			delta.col_z().length() ) < tolerance );
	}

	/// @brief convert a global reference (lab) frame vector to our local (stub) frame
	Vector
	global2local( Vector const & xyz ) const
	{
		return M.transposed() * ( xyz - v );
	}

	/// @brief convert a local reference (stub) frame vector to the global (lab) frame
	Vector
	local2global( Vector const & xyz ) const
	{
		return M * xyz + v;
	}

	/// @brief  Build a vector in the global lab frame from the spherical coords used in the atomtree
	///
	/// @details theta is the angle between the postive x and the vector (0<=theta<=pi),
	///phi is the angle between the y-z plane projection of the vector and the positive y (0<=theta<=2*pi),
	///d is the length of the vector
	///
	/// @note  These are non-standard in the choice of axis for historical reasons --PB
	Vector
	spherical( Real const phi, Real const theta, Real const d ) const
	{
		Real const d_sin_theta( d * std::sin( theta ) );
		return local2global( Vector(
			d * std::cos( theta ),
			d_sin_theta * std::cos( phi ),
			d_sin_theta * std::sin( phi ) ) );
	}

	/// @brief  Build stubatom coords that would yield this stub
	Vector
	build_fake_xyz( Size const index ) const;

	/// @brief  Build stubatom coords that would yield this stub
	Vector
	build_fake_xyz( Size const index, Real const length1, Real const length2, Real const angle_degrees ) const;

	bool operator == ( Stub const & rhs ) const;

public: // Fields
	/// @brief coord frame by 3x3 matrix, each column is a unit vector
	Matrix M;

	/// @brief Coordinate frame rotation matrix.
	Matrix rotation() const { return M; }

	/// @brief center point by a vector
	Vector v;

	/// @brief Coordinate frame center.
	Vector center() const { return v; }
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // Stub

/// @brief root squared deviation between two stubs
inline
double
distance(
	Stub const & a,
	Stub const & b
)
{
	using namespace numeric;
	return std::sqrt(
		a.M.col_x().distance_squared(b.M.col_x() ) +
		a.M.col_y().distance_squared(b.M.col_y() ) +
		a.M.col_z().distance_squared(b.M.col_z() ) +
		a.v.distance_squared(        b.v         )
	);
}


std::ostream &
operator<<( std::ostream & os, Stub const & a );

/// @brief Global default stub
extern Stub default_stub;


} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_Stub_HH
