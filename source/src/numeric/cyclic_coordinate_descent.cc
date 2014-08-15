// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/cyclic_coordinate_descent.cc
/// @brief Compute the angle that minimizes the deviation between a set of atoms
/// @author Brian D. Weitzner
/// @author Labonte <JWLabonte@jhu.edu>
/// @author Phil Bradley

#include <numeric/cyclic_coordinate_descent.hh>
#include <numeric/conversions.hh>

namespace numeric {

typedef xyzVector< Real > Vector;
typedef Real Length;


/// @param <F>: the coordinates of the fixed target atoms
/// @param <M>: the coordinates of the moving positions to be overlapped with the target atoms
/// @param <theta_hat>: axis vector of the torsion angle
/// @param <alpha>: empty angle to be calculated
/// @param <S>: empty deviation to be calculated
///
/// @details The objective of an individual cyclic coordinate descent (CCD) move is to minimize the deviation between
/// a set of points that should perfectly superimpose. The deviation squared (S) can be expressed as:
///
///     S = Sum(r^2 + f^2) - 2 Sum[r(f_vector dot r_hat)] cos theta - 2 Sum[r(f_vector dot s_hat)] sin theta
///
/// The derivative of S with respect to theta (the angle about the rotation axis):
///
///     dS/dtheta = 2 Sum[r(f_vector dot r_hat)] sin theta - 2 Sum[r(f_vector dot s_hat)] cos theta
///
/// Setting dS/dtheta to zero gives the minimal value of theta, which we call alpha:
///
///     tan alpha = Sum[r(f_vector dot s_hat] / Sum[r(f_vector dot r_hat]
///
/// If we define...
///     a = Sum(r^2 + f^2)
///     b = 2 Sum[r(f_vector dot r_hat)]
///     c = 2 Sum[r(f_vector dot s_hat)]
/// then S can be rewritten:
///     S = a - b cos alpha - c sin alpha
/// and we can express alpha as
///     tan alpha = c / b
void ccd_angle(
	utility::vector1< xyzVector< Real > > const & F,
	utility::vector1< xyzVector< Real > > const & M,
	xyzVector< Real > const & axis_atom,
	xyzVector< Real > const & theta_hat,
	Real & alpha,
	Real & S)
{
	// Make sure theta_hat is a unit vector
	// Should we try to recover in the case where it is not?
	assert( theta_hat.is_unit( 1e-3 ) );

	// Accumulate the components of the CCD equation.  These values are explained below.
	Real aa( 0.0 );
	Real bb( 0.0 );
	Real cc( 0.0 );

	Size n_overlap_positions( min( F.size(), M.size() ) );

	// For each position to overlap, calculate the variables required for computing the best angle, alpha.
	for ( Size i = 1; i <= n_overlap_positions; ++i ) {
		// O(i) is the origin of the vector from the axis to coordinates F(i) and M(i).
		Vector const O( axis_atom + dot( M[ i ] - axis_atom, theta_hat ) * theta_hat );

		// f(i) is the vector from O(i) to F(i).
		Vector const f_vector( F[ i ] - O );
		Length const f_length( f_vector.length() );

		// r(i) is the vector from O(i) to M(i).
		Vector const r_vector( M[ i ] - O );
		Length const r_length( r_vector.length() );

		aa += ( r_length * r_length ) + ( f_length * f_length );

		// If r_length is zero, M(i) is O(i) and its position cannot change; move to the next i.
		if ( r_length < 1e-6 ) { continue; }  // TODO: Find & add note about PDB precision.
		Vector const r_hat = r_vector.normalized();  // unit vector for r

		// s_hat(i) is the unit vector normal to the axis and r_hat(i); orientation choice comes here!
		Vector const s_hat( cross( theta_hat, r_hat ) );

		// aa += r_length * r_length + f_length * f_length;

		bb += 2 * r_length * dot( f_vector, r_hat );

		cc += 2 * r_length * dot( f_vector, s_hat );

		// some checks on orthonormality
		assert( abs( dot( theta_hat, r_hat ) ) < 1e-3 );
		assert( abs( dot( theta_hat, s_hat ) ) < 1e-3 );
		assert( abs( dot( s_hat, r_hat ) ) < 1e-3 );
		assert( s_hat.is_unit( 1e-3 ) );
		assert( r_hat.is_unit( 1e-3 ) );

		assert( abs( dot( M[ i ], theta_hat ) - dot( O, theta_hat ) ) < 1e-3 );
		assert( abs( dot( r_vector, theta_hat ) ) < 1e-3 );
	}

	// Derivation shown above.
	alpha = atan2( cc, bb );

	// Derivation shown above.
	S = aa - bb * cos( alpha ) - cc * sin( alpha );

	alpha = conversions::degrees( alpha );
}

} // numeric
