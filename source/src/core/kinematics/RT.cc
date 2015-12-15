// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/RT.cc
/// @brief  Rotation + translation class
/// @author Phil Bradley


// Unit headers
#include <core/kinematics/RT.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>

// C++ Headers
#include <iostream>
#include <string>

#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

#endif // SERIALIZATION


namespace core {
namespace kinematics {

///////////////////////////////////////////////////////////////////////////////
/// @details rb are the small changes in translation and rotation, written in up-stream (jump_start)
/// rb_center is the center of rotation, written in the down-stream (jump_end)

// coordinate system NOTE: downstream!!!!!
//
void
RT::fold_in_rb_deltas(
	utility::vector1<Real> const & rb,
	Vector const & rb_center
)
{
	using namespace numeric;

	debug_assert( Size(rb.size()) == 6 );
	// simple: rotation gets multiplied by Rzyx,
	// translation (t) goes to center + Rzyx(t-center) + rb_trans

	// find position of rb-center in jump_start coord sys
	// m2 -> m1: R * V
	Vector new_center(translation + rotation * rb_center);


	// create a transformation matrix from the 3 rb angles
	Matrix const Rzyx(
		z_rotation_matrix_degrees( rb[6] ) * (
		y_rotation_matrix_degrees( rb[5] ) *
		x_rotation_matrix_degrees( rb[4] ) )
	);

	rotation.left_multiply_by( Rzyx );


	Vector rb_trans( rb[1], rb[2], rb[3] ); // load the first 3
	translation = new_center + Rzyx * ( translation - new_center ) + rb_trans;
} // fold_in_rb_deltas


///////////////////////////////////////////////////////////////////////////////
bool
RT::ortho_check()
const
{
	Real const tolerance( 1e-3 );
	Matrix delta( rotation * rotation.transposed() - Matrix::identity() );
	return ( ( delta.col_x().length() +
		delta.col_y().length() +
		delta.col_z().length() ) < tolerance );
}
//  // debug orthogonality
//  Real dev(0.0);
//  for ( int i=1; i<=3; ++i ) {
//   for ( int j=1; j<=3; ++j ) {
//    Real f =
//     rotation(1,i) * rotation(1,j) +
//     rotation(2,i) * rotation(2,j) +
//     rotation(3,i) * rotation(3,j);
//    if ( i==j ) dev += std::abs(1.0-f);
//    else dev += std::abs(f);
//   }
//  }
//  return ( dev < 0.01 );
// }


/////////////////////////////////////////////////////////////////////////////
// @details these two extractors must be inverses: ie if we write an RT
// out with one and then read it back in using the other the
// new RT should be the same
// this is critical for pose silent input/output to work correctly.

std::ostream &
operator <<(
	std::ostream & os,
	const RT & rt
)
{
	debug_assert( rt.ortho_check() );
	os << "RT "; // olange: removed whitespace before RT --> gives problem in silent files.
	for ( int i = 1; i <= 3; ++i ) {
		for ( int j = 1; j <= 3; ++j ) {
			//os << F(9,3,rt.rotation(i,j)); // output by row
			os << rt.rotation(i,j) << ' '; // chu -- more digits for precision
		}
	}
	for ( int i = 1; i <= 3; ++i ) {
		//os << F(9,3,rt.translation(i) );
		os << rt.translation(i) << ' '; // chu -- more digits for precision
	}
	os << ' ';
	return os;
}

///////////////////////////////////////////////////////////////////////////////
std::istream &
operator >>(
	std::istream & is,
	RT & rt
)
{
	std::string tag;
	is >> tag;
	if ( !is.fail() && tag == "RT" ) {
		for ( int i = 1; i <= 3; ++i ) {
			for ( int j = 1; j <= 3; ++j ) {
				is >> rt.rotation(i,j); // input by row
			}
		}
		for ( int i = 1; i <= 3; ++i ) {
			is >> rt.translation(i);
		}
		if ( !is.fail() ) {
			if ( !rt.ortho_check() ) {
				std::cout << "RT failed orthogonality check!" << std::endl;
				is.setstate( std::ios_base::failbit );
			}
		}
	}
	return is;
}


} // namespace kinematics
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::kinematics::RT::save( Archive & arc ) const {
	arc( CEREAL_NVP( rotation ) ); // Matrix
	arc( CEREAL_NVP( translation ) ); // Vector
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::kinematics::RT::load( Archive & arc ) {
	arc( rotation ); // Matrix
	arc( translation ); // Vector
}

SAVE_AND_LOAD_SERIALIZABLE( core::kinematics::RT );
#endif // SERIALIZATION
