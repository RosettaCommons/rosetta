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

// Project headers
#include <core/kinematics/jacobian/util.hh>

// Package headers
#include <basic/Tracer.hh>
#include <numeric/xyz.functions.hh>
#include <core/kinematics/RT.hh>

namespace core {
namespace kinematics {
namespace jacobian {

static basic::Tracer TR( "core.kinematics.jacobian.util" );

void
weigh_columns_inversely_squared(utility::vector1< numeric::MathVector<core::Real> > & input_cols){
	// Input is a 1-vector with MathVectors that all represent the same target quantity.
	// Output is a 1-vector with MathVectors where the contribution to this quantity is spread out over the vectors
	// so that the sum of their contribution is the same target quantity.
	// In this implementation each contribution is weighed as the squared inverse of its size, so that the target quantity is
	// reached with 'minimal effort'.

	// total weight factor
	core::Real total_weight = 0;
	// 1-vector with individual weights
	utility::vector1< core::Real > col_weights( input_cols.size(), 0.0 );

	for ( core::Size i=1; i <= input_cols.size(); ++i ) {
		// calculate norm
		core::Real col_norm = input_cols[i].square_norm();
		// exclude columns that have a zero (or very low) norm
		if ( col_norm > margin_from_zero_in_jacobian_util ) {
			// correct for number of DoFs in case it's < 6 DoF
			col_weights[i] = 1 / col_norm;
		}
		// increment weight
		total_weight += col_weights[i];
	}

	// only update cols if total weight was large enough
	if ( total_weight > margin_from_zero_in_jacobian_util ) {
		for ( core::Size i = 1; i <= input_cols.size(); ++i ) {
			input_cols[i] *= col_weights[i] / total_weight;
		}
	} // otherwise keep the input_cols as they were
}

void
weigh_columns_by_energy_grad(utility::vector1< numeric::MathVector<core::Real> > & input_cols, utility::vector1< numeric::MathVector<core::Real> > const & gradient_cols, core::Real const scaling){
	// First input is a 1-vector with MathVectors of differential torsion angles, which all represent the same differential
	// dual Cartesian vector (screw).
	// Second input is a 1-vector with MathVectors that represent the energy gradient for each represented DoF.
	// Output is a 1-vector with MathVectors where the contribution to this quantity is spread out over the vectors
	// so that the sum of their contribution is the same target quantity.
	// In this implementation each contribution is weighed according to the -dE/norm(col), capped at a minimum value
	// of 0.01, so that positive gradients get a low weighing.

	// total weight factor
	core::Real total_weight = 0;
	// 1-vector with individual weights
	utility::vector1< core::Real > col_weights( input_cols.size(), 0.0 );

	for ( core::Size i=1; i <= input_cols.size(); ++i ) {
		// calculate norm of vars change
		core::Real col_norm = input_cols[i].norm();
		// only include columns that have a non-zero norm
		if ( col_norm > margin_from_zero_in_jacobian_util ) { // inverted sigmoid function with scaling factor
			// calculate net energy descent for this module
			core::Real dE(0.0);
			for ( core::Size j=0; j < 6; ++j ) {
				dE += input_cols[i](j) * gradient_cols[i](j);
			}
			col_weights[i] = 1.0 / (1.0 + std::exp(dE / col_norm * scaling) );
		}
		// increment weight
		total_weight += col_weights[i];
	}

	// only update cols if total weight was large enough
	if ( total_weight > margin_from_zero_in_jacobian_util ) {
		for ( core::Size i = 1; i <= input_cols.size(); ++i ) {
			input_cols[i] *= col_weights[i] / total_weight;
		}
	} // otherwise keep the input_cols as they were
}

numeric::xyzVector<core::Real>
project_ZYX_euler_on_cart(numeric::xyzVector<core::Real> const euler_angles, numeric::xyzVector<core::Real> const euler_error) {
	numeric::xyzVector<core::Real> Cart_rot_error =
		// projection of error vector about instantaneous roll axis
		numeric::z_rotation_matrix_radians(euler_angles(1)) * numeric::y_rotation_matrix_radians(euler_angles(2)) *
		numeric::xyzVector<core::Real>(euler_error(3), 0, 0) +
		// plus projection of error vector about instantaneous pitch axis
		numeric::z_rotation_matrix_radians(euler_angles(1)) *
		numeric::xyzVector<core::Real>(0, euler_error(2), 0) +
		// plus projection of error vector about instantaneous yaw axis
		numeric::xyzVector<core::Real>(0, 0, euler_error(1));

	return Cart_rot_error;
}


/// @brief Extract the error twist (linearized rotations) from the current and target position and orientation (in Euler angles)
ModuleType1::Screw
calculate_error_twist( core::kinematics::RT const RT_current, core::kinematics::RT const RT_target ) {
	// create RT object of the error
	core::kinematics::RT const RT_error(RT_current, RT_target);
	// create HomogeneousTransform so can use euler projection
	numeric::HomogeneousTransform<core::Real> const H_error(RT_error.get_rotation(), RT_error.get_translation());
	// obtain ZYX euler angles of error matrix
	numeric::xyzVector<core::Real> const euler_error = H_error.euler_angles_ZYX_rad();
	// project euler error on Cartesian. This is technically flawed, and in theory only works for small angles
	numeric::xyzVector<core::Real> const Cart_rot_error_local = project_ZYX_euler_on_cart(euler_error, euler_error);
	// Transform linearized Euler angle to global reference frame
	numeric::xyzVector<core::Real> const Cart_rot_error(RT_current.get_rotation() * Cart_rot_error_local);
	// get linear error by transforming the linear error into the reference frame connected to the ref atom
	numeric::xyzVector<core::Real> const Cart_lin_error = RT_target.get_translation() - RT_current.get_translation() -
		Cart_rot_error.cross_product(RT_current.get_translation());
	// put results in screw error
	core::kinematics::jacobian::ModuleType1::Screw dT;
	dT.first = Cart_rot_error;
	dT.second = Cart_lin_error;

	// return answer
	return dT;
}

} //jacobian
} //kinematics
} //core
