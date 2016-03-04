// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/RT.hh
/// @brief  Rotation + translation class
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_RT_hh
#define INCLUDED_core_kinematics_RT_hh


// Unit headers
#include <core/kinematics/RT.fwd.hh>

// Package headers
#include <core/kinematics/Stub.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// Utility headers
#include <utility/vector1.hh>

//#include <utility/pointer/ReferenceCount.hh>


namespace core {
namespace kinematics {


/// @brief Rotation + translation class
/// @note  Used in jumps
//class RT : public utility::pointer::ReferenceCount {
class RT {

public: // Types

	typedef  numeric::xyzVector< Real >  Vector; // DOUBLE!
	typedef  numeric::xyzMatrix< Real >  Matrix; // DOUBLE!

	// Types to prevent compile failure when std::distance is in scope
	typedef  void  iterator_category;
	typedef  void  difference_type;

public:

	/// @brief Default constructor
	inline
	RT():
		rotation( Matrix::identity() ),
		translation( 0.0 )
	{}

	/// @brief matrix/vector constructor
	inline
	RT( Matrix const & rot, Vector const & vec ):
		rotation( rot ),
		translation( vec )
	{}

	/// @brief construct from two stubs
	RT( Stub const & stub1, Stub const & stub2 )
	{
		this->from_stubs( stub1, stub2 );
	}

	/// @brief copy constructor
	RT ( RT const & src ):
		//ReferenceCount(),
		rotation( src.rotation ),
		translation( src.translation )
	{}

	/// @brief update from stubs
	void
	from_stubs(
		Stub const & stub1,
		Stub const & stub2
	);

	/// @brief reverse the "view"
	void
	reverse()
	{
		// new rotation is the inversed(transposed) old rotation
		rotation.transpose();

		// new translation is the negated old tranlation in m1 frame rewritten
		// in m2 frame (left multiplied by old_rotation^T = new_rotation)
		translation = rotation * translation.negate();
	}

	/// @brief reverse the "view"
	void
	inverse()
	{
		reverse();
	}

	/// @brief return to default-constructed state
	void
	reset()
	{
		translation.zero();
		rotation.to_identity();
	}

	/// @brief return to default-constructed state (another name)
	void
	identity_transform()
	{
		translation.zero();
		rotation.to_identity();
	}

	/// @brief set the tranlsation
	void
	set_translation( Vector const & t )
	{
		translation = t;
	}

	/// @brief set the rotation
	void
	set_rotation( Matrix const & r )
	{
		rotation = r;
	}

	/// @brief get the rotation
	Matrix const &
	get_rotation() const
	{
		return rotation;
	}

	/// @brief get the translation
	Vector const &
	get_translation() const
	{
		return translation;
	}

	Real
	distance_squared( RT const & b ) const
	{
		return ( translation    .distance_squared( b.translation ) +
			rotation.col(1).distance_squared( b.rotation.col(1) ) +
			rotation.col(2).distance_squared( b.rotation.col(2) ) +
			rotation.col(3).distance_squared( b.rotation.col(3) ) );
	}


	/// @brief update the transform to include small additional rigid-body rotations and translations
	/// @note PHIL: IT WOULD BE GOOD TO ELIMINATE ARGUMENT ARRAYS IN MINI BY USE OF
	///  APPROPRIATE LAYERED DATA STRUCTURES
	void
	fold_in_rb_deltas(
		utility::vector1<Real> const &  rb,
		Vector const & rb_center
	);


	/// @brief should be inverse of RT::from_stubs
	void
	make_jump(
		Stub const & stub1,
		Stub & stub2 ) const;

	/// @brief output operator
	friend std::ostream & operator <<( std::ostream & os, RT const & rt );
	/// @brief input operator
	friend std::istream & operator >>( std::istream & is, RT & rt );

	/// @brief copy operator
	RT & operator =( RT const & src ) {
		rotation = src.rotation;
		translation = src.translation;
		return *this;
	}

	/// @brief check for orthogonality
	bool ortho_check() const;

	friend
	inline
	Real
	distance( RT const & a, RT const & b );

private: // Fields
	/// rotation matrix, written in stub1 frame
	Matrix rotation; // 3x3
	/// tranlsation vector, written in stub1 frame
	Vector translation; // 3

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // RT


///////////////////////////////////////////////////////////////////////////////
/// @brief root squared devitation of two RTs
inline
Real
distance( RT const & a, RT const & b )
{
	using namespace numeric;

	return
		std::sqrt( a.rotation.col(1).distance_squared( b.rotation.col(1) ) +
		a.rotation.col(2).distance_squared( b.rotation.col(2) ) +
		a.rotation.col(3).distance_squared( b.rotation.col(3) ) +
		a.translation    .distance_squared( b.translation ) );
}

} // kinematics
} // core


#endif // INCLUDED_core_kinematics_RT_HH
