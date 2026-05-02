// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/crick_equations/BundleParams_derivatives.hh
/// @brief  Analytical derivatives of the Crick equations with respect to bundle parameters.
/// @details Computes the Jacobian dXYZ/dParam for each perturbable Crick parameter,
/// enabling gradient-based minimization in parametric space.
/// @author Andy Watkins

#ifndef INCLUDED_numeric_crick_equations_BundleParams_derivatives_hh
#define INCLUDED_numeric_crick_equations_BundleParams_derivatives_hh

#include <numeric/types.hh>
#include <numeric/xyzVector.hh>

namespace numeric {
namespace crick_equations {

/// @brief Derivatives of Crick XYZ coordinates with respect to each perturbable parameter.
struct CrickDerivatives {
	xyzVector<Real> dXYZ_dr0;
	xyzVector<Real> dXYZ_domega0;
	xyzVector<Real> dXYZ_ddelta_omega0;
	xyzVector<Real> dXYZ_ddelta_omega1;
	xyzVector<Real> dXYZ_ddelta_t;
};

/// @brief Compute analytical derivatives of the Crick bundle equations.
/// @details Currently supports the epsilon=1 (circular cross-section) case only.
/// @param[in] t Position along the major helix (residue index, can be fractional)
/// @param[in] r0 Major helix radius (Angstroms)
/// @param[in] omega0 Major helix angular frequency (radians per residue)
/// @param[in] delta_omega0 Major helix phase offset (radians)
/// @param[in] r1 Minor helix radius for this atom (Angstroms)
/// @param[in] omega1 Minor helix angular frequency (radians per residue)
/// @param[in] z1 Minor helix rise per residue (Angstroms)
/// @param[in] delta_omega1 Minor helix phase offset for this atom (radians)
/// @param[in] delta_z1 Axial offset for this atom (Angstroms)
/// @param[out] failed Set to true if the derivative computation fails
/// @return CrickDerivatives struct with dXYZ/dParam for each parameter
CrickDerivatives compute_crick_derivatives(
	Real t,
	Real r0,
	Real omega0,
	Real delta_omega0,
	Real r1,
	Real omega1,
	Real z1,
	Real delta_omega1,
	Real delta_z1,
	bool & failed
);

} // namespace crick_equations
} // namespace numeric

#endif
