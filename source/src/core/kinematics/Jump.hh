// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/Jump.hh
/// @brief  Kinematics Jump class
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_Jump_hh
#define INCLUDED_core_kinematics_Jump_hh


// Unit headers
#include <core/kinematics/Jump.fwd.hh>

// Package Headers
#include <core/kinematics/RT.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

// Utility Headers
#include <utility/vector1.hh>


namespace core {
namespace kinematics {

// This needs to be a function to avoid the static intitilization order fiasco
utility::vector1<Real> const & ZERO();

/// @brief an object which makes rigid-body transformation with translational and rotational perturbation
///
/// See @ref atomtree_overview "AtomTree overview and concepts" for details.
///
class Jump {
public: // Types
	static const int TRANS_X = 1;
	static const int TRANS_Y = 2;
	static const int TRANS_Z = 3;
	static const int ROT_X = 4;
	static const int ROT_Y = 5;
	static const int ROT_Z = 6;

	typedef  numeric::xyzVector< Real >  Vector; // DOUBLE!
	typedef  numeric::xyzMatrix< Real >  Matrix; // DOUBLE!

public:
	/// @brief construction
	Jump ():
		rt_(),
		rb_delta( 2, ZERO() ),
		rb_center( 2, Vector(0.0) ),
		invert_upstream_(false),
		invert_downstream_(false)
	{}

	/// @brief constructor with only RT
	Jump( const RT & src_rt ):
		rt_( src_rt ),
		rb_delta ( 2, ZERO() ),
		rb_center( 2, Vector(0.0) ),
		invert_upstream_(false),
		invert_downstream_(false)
	{}

	/// @brief get RT from two stubs and ZERO rb_delta
	Jump ( Stub const & stub1, Stub const & stub2 ):
		rt_( stub1, stub2 ),
		rb_delta( 6, ZERO() ),
		rb_center( 2, Vector(0.0) ),
		invert_upstream_(false),
		invert_downstream_(false)
	{}

	/// @brief copy constructor
	Jump ( const Jump & )= default;

	/// @brief copy operator
	Jump &
	operator =( Jump const & ) = default;

	/// @brief get a jump from two stubs
	void
	from_stubs(
		Stub const & stub1,
		Stub const & stub2
	);

	/// @brief get a jump from a bond cst definition
	void
	from_bond_cst(
		utility::vector1< Vector > & atoms,
		utility::vector1< Real > const & csts
	);

	/// @brief make a jump from stub1 to stub2
	void
	make_jump(
		Stub const & stub1,
		Stub & stub2
	) const;

	/// @brief check RT's orthogonality
	bool ortho_check() const;

	/// @brief check whether there is rb_delta not being transformed.
	bool nonzero_deltas() const;

	/// @brief reset RT, rb_delta and rb_center
	void reset();

	/// @brief translate along a randomly chosen vector by dist_in
	void random_trans( const float dist_in );

	/// @brief Given the desired magnitude of the translation and rotation
	/// components, applies Gaussian perturbation to this jump.
	/// Return the move that was applied
	utility::vector1<Real> gaussian_move(int const dir, float const trans_mag, float const rot_mag);

	// @brief make a gaussian move with one selected rb dof
	void
	gaussian_move_single_rb(
		int const dir,
		float const mag,
		int rb
	);

	/// @brief make a rotation "matrix" about the center "center"
	void
	rotation_by_matrix(
		Stub const & stub,
		Vector const & center, //in xyz frame
		Matrix const & matrix //in xyz frame
	);

	/// @brief make a rotation of alpha degrees around axix and center
	void
	rotation_by_axis(
		Stub const & stub,
		Vector const & axis,
		Vector const & center, //in xyz frame
		float const alpha
	);

	/// @brief make a translation along axis by dist
	void
	translation_along_axis(
		Stub const & stub,
		Vector const & axis,
		float const dist
	);

	/// @brief reset to identity matrix, 0 translation
	void
	identity_transform();

	/// @brief change the direction of a jump, e.g, upstream stub becomes downstream stub now.
	void
	reverse();

	/// @brief Return a new jump that is the inverse transformation
	Jump
	reversed() const;

	/// @brief get rotation matrix
	inline
	Matrix const &
	get_rotation() const;

	/// @brief get translation vector
	inline
	Vector const &
	get_translation() const;

	/// @brief set rotation matrix
	void
	set_rotation( Matrix const & R_in );

	/// @brief set translation vector
	void
	set_translation( Vector const & t );

	/// @brief transform small changes in translation and rotation into jump
	void
	fold_in_rb_deltas(); // call after we're done minimizing

	/// @brief return the RT modeled by this jump --> makes new Jump, calls fold_in_rb_delte and returns RT
	RT rt() const;

	/// @brief get rb_delta by direction
	inline
	utility::vector1<Real> const &
	get_rb_delta( int const dir ) const;

	/// @brief get rb_delta by direction and rb_number
	inline
	Real
	get_rb_delta( int const rb_no, int const dir ) const;

	/// @brief set rb_delta by direction and rb_number
	void
	set_rb_delta( int const rb_no, int const dir, Real const value );

	/// @brief set rb_deltas by direction
	void
	set_rb_deltas( int const dir, utility::vector1<Real> const & );

	/// @brief get rb_center by direction
	inline
	Vector const &
	get_rb_center( int const dir ) const;

	/// @brief set rb_center by direction
	void
	set_rb_center(
		const int dir,
		Stub const & stub,
		Vector const & center
	);

	inline bool
	get_invert_upstream( ) const { return invert_upstream_; }
	inline bool
	get_invert_downstream( ) const { return invert_downstream_; }

	// set the upstream or downstream inversion
	// note that these need the upstream stub in order to properly compute the
	//   new rotation matrix
	void
	set_invert( bool upstream, bool downstream);

	bool operator==( Jump const& ) const;

	bool operator!=( Jump const& other ) const { return !operator==( other ); }

	/// @brief stream output operator
	friend std::ostream & operator <<( std::ostream & os, const Jump & jump );
	/// @brief stream input  operator
	friend std::istream & operator >>( std::istream & is, Jump & jump );
	/// @brief RT root squared deviation
	friend Real distance( Jump const & a_in, Jump const & b_in );

	// convenience
	static Matrix mirror_z_transform;

private:
	// private methods


	/// @brief get the lookup index into rb_delta and rb_center from the direction dir.
	/** dir should be 1 or -1 */

	inline
	int
	rb_index( int const dir ) const
	{
		debug_assert( dir == 1 || dir == -1 );
		return ( dir == 1 ? 1 : 2 );
	}


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	/// @brief translation and rotation
	///
	RT rt_;

	/// @brief changes to translation and rotation
	/// @details 6x2 table, for each of the two folding directions, the first three are for translations
	/// along xyz, and the next three are for rotatation around xyz axes.
	utility::vector1< utility::vector1<Real> > rb_delta; // 6x2

	/// @brief the center around which the rotation is performed
	/// @details 3x2 table, for each of the two folding directions, the rotation center is written in the
	/// local frame of the downstream stub.
	utility::vector1< Vector > rb_center; // 3x2

	/// @brief are the upstream or downstream residues an inverted coordinate system (Z <- -X cross Y)
	bool invert_upstream_, invert_downstream_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // Jump


///////////////////////////////////////////////////////////////////////////////
/// inline definitions for Jump class

inline
Jump::Matrix const &
Jump::get_rotation() const
{
	return rt_.get_rotation();
}

///////////////////////////////////////////////////////////////////////////////
inline
Jump::Vector const &
Jump::get_translation() const
{
	return rt_.get_translation();
}

///////////////////////////////////////////////////////////////////////////////
// this one is private

///////////////////////////////////////////////////////////////////////////////
//
inline
const utility::vector1<Real> &
Jump::get_rb_delta( int const dir ) const
{
	return rb_delta[ rb_index( dir ) ];
}

///////////////////////////////////////////////////////////////////////////////
//
inline
Real
Jump::get_rb_delta( int const rb_no, int const dir ) const
{
	debug_assert( dir == 1 || dir == -1 );
	return rb_delta[ rb_index( dir ) ][ rb_no ] ;
}

///////////////////////////////////////////////////////////////////////////////
//
inline
Jump::Vector const &
Jump::get_rb_center( int const dir ) const
{
	debug_assert( dir == 1 || dir == -1 );
	return rb_center[ rb_index( dir ) ];
}

/// @brief compare the difference of two jumps in term of the translation (dist) and rotational angle(theta)
void
jump_distance(
	Jump const & a_in,
	Jump const & b_in,
	Real & dist,
	Real & theta
);


/// @brief RT root squared deviation
Real distance( Jump const & a_in, Jump const & b_in );


} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_Jump_HH
