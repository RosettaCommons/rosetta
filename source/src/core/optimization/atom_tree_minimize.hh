// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/atom_tree_minimize.hh
/// @brief  Atom tree minimization functions
/// @author Phil Bradley

#ifndef INCLUDED_core_optimization_atom_tree_minimize_hh
#define INCLUDED_core_optimization_atom_tree_minimize_hh


// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/NumericalDerivCheckResult.fwd.hh>
#include <core/optimization/DOF_Node.fwd.hh>
#include <core/kinematics/tree/Atom.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// ObjexxFCL headers

#include <core/optimization/MinimizerMap.fwd.hh>
#include <core/optimization/Multifunc.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace optimization {


void
atom_tree_dfunc(
	pose::Pose & pose,
	MinimizerMap & min_map,
	scoring::ScoreFunction const & scorefxn,
	Multivec const & vars,
	Multivec & dE_dvars
);


void
atom_tree_get_atompairE_deriv(
	pose::Pose & pose,
	MinimizerMap & min_map,
	scoring::ScoreFunction const & scorefxn
);

Real
torsional_derivative_from_cartesian_derivatives(
	kinematics::tree::Atom const & atom,
	optimization::DOF_Node const & dof_node,
	Real dof_deriv, // derivatives applied directly to this DOF, if any
	Real torsion_scale_factor
);

/// @brief Numeric deriv check for Multifuncs other than the AtomTreeMultifunc.
SimpleDerivCheckResult
simple_numeric_deriv_check(
	Multifunc const & func,
	Multivec const & start_vars,
	Multivec const & dE_dvars,
	bool send_to_stdout,
	bool verbose,
	Size nsteps = 5
);

void
numerical_derivative_check(
	MinimizerMap const & min_map,
	Multifunc const & func,
	Multivec const & start_vars,
	Multivec const & dE_dvars,
	NumericalDerivCheckResultOP deriv_check_result,
	bool const verbose // = true
);


// Real
// calculate_direct_dof_derivatives(
//  DOF_Node const & tor,
//  pose::Pose const & pose,
//  scoring::ScoreFunction const & scorefxn,
//  ObjexxFCL::FArray2D< Real > const & dunbrack_deriv // currently this is pre-computed
// );


} // namespace optimization
} // namespace core

#endif // INCLUDED_core_optimization_atom_tree_minimize_HH
