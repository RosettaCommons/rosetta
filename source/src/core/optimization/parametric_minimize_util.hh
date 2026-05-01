// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/parametric_minimize_util.hh
/// @brief  Utility functions for minimizing over parametric DOFs (e.g. Crick parameters for helical bundles).
/// @author Andy Watkins (andy.watkins2@gmail.com)

#ifndef INCLUDED_core_optimization_parametric_minimize_util_hh
#define INCLUDED_core_optimization_parametric_minimize_util_hh

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <set>
#include <string>

namespace core {
namespace optimization {

/// @brief Information about a single parametric DOF to be minimized.
struct ParametricDOFInfo {
	Size params_set_index;    // Index into pose.conformation().parameters_set()
	Size params_index;        // Index into ParametersSet::parameters()
	Size param_enum;          // Parameter enum value (e.g., BPC_r0, BBPC_r0)
	Size helix_start_resid;   // First residue of this helix/strand in the pose
	Size helix_end_resid;     // Last residue of this helix/strand in the pose
	bool is_global;           // True if this param is global for the ParametersSet
	Real scale_factor;        // Scale factor for DOF vector (default 1.0)
	std::string param_name;   // Parameter name (e.g., "r0", "omega0") for Jacobian dispatch

	ParametricDOFInfo() :
		params_set_index( 0 ),
		params_index( 0 ),
		param_enum( 0 ),
		helix_start_resid( 0 ),
		helix_end_resid( 0 ),
		is_global( false ),
		scale_factor( 1.0 ),
		param_name( "" )
	{}
};

/// @brief Enumerate all parametric DOFs from a pose that has ParametersSets.
/// @details For each ParametersSet, for each Parameters, for each perturbable
/// RealValuedParameter (those with can_be_perturbed()==true), creates a ParametricDOFInfo.
/// For parameters marked as global_for_parameters_set, only one DOF is created
/// (associated with the first Parameters object).
void enumerate_parametric_dofs(
	pose::Pose const & pose,
	utility::vector1< ParametricDOFInfo > & dof_infos
);

/// @brief Get the set of residue indices that are under parametric control.
/// @details Reads the residue lists from all Parameters objects in all ParametersSets.
std::set< Size > get_parametric_residues(
	pose::Pose const & pose
);

/// @brief Get the current value of a parametric DOF from the pose.
Real get_parametric_dof_value(
	pose::Pose const & pose,
	ParametricDOFInfo const & info
);

/// @brief Set a parametric DOF value in the pose's ParametersSet.
void set_parametric_dof_value(
	pose::Pose & pose,
	ParametricDOFInfo const & info,
	Real value
);

/// @brief Rebuild backbone coordinates for all parametric elements in the pose.
/// @details Iterates over all ParametersSets, reconstructing backbone atom positions
/// from current parameter values via the ParametrizationCalculator.
/// This is a potentially expensive operation (full Crick equation evaluation).
void rebuild_parametric_backbone(
	pose::Pose & pose
);

} // namespace optimization
} // namespace core

#endif // INCLUDED_core_optimization_parametric_minimize_util_hh
