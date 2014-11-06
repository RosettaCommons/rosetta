// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file BundleParams.hh.
/// @brief Functions implementing the Crick equations for a straight helix (NOT a helical bundle).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_numeric_crick_equations_BundleParams_hh
#define INCLUDED_numeric_crick_equations_BundleParams_hh

//C++ headers
#include <math.h>

// Unit headers
#include <numeric/types.hh>
#include <numeric/xyzVector.hh>

namespace numeric {
	namespace crick_equations {

	/// @brief Calculates P, the repeat distance of the major helix in the z-direction.
	/// @details Sets "failed" to TRUE if and only if the value of P could not be calculated.
	Real PP (
		Real const &r0,
		Real const &omega0,
		Real const &z1,
		bool &failed
	);

	/// @brief Calculates P*omega0 (the repeat distance of the major helix in the z-direction times the major helix
	/// turn per residue).
	/// @details Sets "failed" to TRUE if and only if the value of P*omega0 could not be calculated.
	Real P_omega0 (
		Real const &r0,
		Real const &omega0,
		Real const &z1,
		bool &failed
	);

	/// @brief Calculates alpha, the tilt angle of the minor helix.
	///
	Real ALPHA (
		Real const &r0,
		Real const &omega0,
		Real const &P_omega0
	);

	/// @brief Calculates the norm of the gradient vector
	/// @details Needed for adding the small delta_z1 offsets, which are in the direction
	/// of (dx/dt, dy/dt, dz/dt) of the major helix.
	Real gradnorm(
		Real const &r0,
		Real const &omega0,
		Real const &Pomega0,
		Real const &s0,
		Real const &c0
	);

	/// @brief Calculates cos(omega*t+delta_omega)
	///
	Real COSFXN( Real const &t, Real const &omega, Real const&delta_omega );

	/// @brief Calculates sin(omega*t+delta_omega)
	///
	Real SINFXN( Real const &t, Real const &omega, Real const&delta_omega );

	/// @brief Calculates the x-position on the helix of helices
	/// given the Crick parameters for the bundle.
	/// @details Returns failed=true if calculation fails, failed=false otherwise.
	Real X_BUNDLE (
		Real const &t,
		Real const &r0,
		Real const &omega0,
		Real const &delta_omega0,
		//Real const &P,
		Real const &r1,
		Real const &omega1,
		Real const &z1,
		Real const &delta_omega1,
		Real const &delta_z1,
		bool &failed
	);

	/// @brief Calculates the x-position on the helix of helices
	/// given the Crick parameters for the bundle.
	/// @details Returns failed=true if calculation fails, failed=false otherwise.
	Real Y_BUNDLE (
		Real const &t,
		Real const &r0,
		Real const &omega0,
		Real const &delta_omega0,
		//Real const &P,
		Real const &r1,
		Real const &omega1,
		Real const &z1,
		Real const &delta_omega1,
		Real const &delta_z1,
		bool &failed
	);

	/// @brief Calculates the x-position on the helix of helices
	/// given the Crick parameters for the bundle.
	/// @details Returns failed=true if calculation fails, failed=false otherwise.
	Real Z_BUNDLE (
		Real const &t,
		Real const &r0,
		Real const &omega0,
		Real const &delta_omega0,
		//Real const &P,
		Real const &r1,
		Real const &omega1,
		Real const &z1,
		Real const &delta_omega1,
		Real const &delta_z1,
		bool &failed
	);

	/// @brief Calculate the x,y,z coordinates of a point on the helix of
	/// helices given the Crick parameters for the bundle.
	/// @details Not quite as efficient as it could be, but it probably
	/// doesn't matter.  (If I really wanted to optimize this, I'd make sure
	/// that c0, c1, s0, s1, and alpha were all calculated once rather than
	/// thrice.)  Returns failed=true if calculation fails, failed=false otherwise.
	xyzVector <Real> XYZ_BUNDLE (
		Real const &t,
		Real const &r0,
		Real const &omega0,
		Real const &delta_omega0,
		//Real const &P,
		Real const &r1,
		Real const &omega1,
		Real const &z1,
		Real const &delta_omega1,
		Real const &delta_z1,
		bool &failed
	);

	} //namespace crick_equations
} //namespace numeric

#endif
