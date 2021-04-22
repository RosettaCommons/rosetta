// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/kinematics/jacobian/util.hh
/// @brief utility functions that can be applied to (parts of) Jacobian objects or related vectors
/// @author teunhoevenaars (teunhoevenaars@gmail.com)


#ifndef INCLUDED_core_kinematics_jacobian_jacobianoperations_hh
#define INCLUDED_core_kinematics_jacobian_jacobianoperations_hh

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/MathVector.hh>
#include <numeric/MathMatrix.hh>
#include <core/kinematics/jacobian/ModuleType1.hh>
#include <core/kinematics/RT.fwd.hh>

// Class headers
#include <core/types.hh>

namespace core {
namespace kinematics {
namespace jacobian {

/// @brief method that gives weights to different columns (that all represent the same error twist in their respective torsion space),
/// inversely proportional to the norm of the required torsion changes
void
weigh_columns_inversely_squared(utility::vector1< numeric::MathVector<core::Real> > & input_cols);

/// @brief method that gives weights to different columns (that all represent the same error twist in their respective torsion space),
/// according to a sigmoid function that gives greatest weight to input_cols alligned with negative gradient_cols.
/// The shape of the sigmoid can be controlled using the variable 'scaling'
void
weigh_columns_by_energy_grad(utility::vector1< numeric::MathVector<core::Real> > & input_cols, utility::vector1< numeric::MathVector<core::Real> > const & gradient_cols, core::Real const scaling = 0.1);

/// @brief projects intrinsic Euler angles onto extrinsic Euler angles in order XYZ (i.e. projected onto Cartesian unit vectors)
numeric::xyzVector<core::Real>
project_ZYX_euler_on_cart(numeric::xyzVector<core::Real> const euler_angles, numeric::xyzVector<core::Real> const euler_error);

/// @brief calculates the error twist between two RT objects, i.e. the linearised dual vector between them
ModuleType1::Screw
calculate_error_twist(core::kinematics::RT const RT_current, core::kinematics::RT const RT_target);

/// @brief
constexpr core::Real margin_from_zero_in_jacobian_util = 1e-6;

} //jacobian
} //kinematics
} //core

#endif //INCLUDED_core_kinematics_jacobian_jacobianoperations_hh
