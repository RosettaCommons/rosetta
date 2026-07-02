// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/parametric_minimize_util.cc
/// @brief  Utility functions for minimizing over parametric DOFs (e.g. Crick parameters for helical bundles).
/// @author Andy Watkins (andy.watkins2@gmail.com)

// Unit headers
#include <core/optimization/parametric_minimize_util.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <core/conformation/parametric/Parameter.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <core/conformation/Residue.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/memory.hh>

// C++ headers
#include <set>

static basic::Tracer TR( "core.optimization.parametric_minimize_util" );

namespace core {
namespace optimization {

/// @brief Enumerate all parametric DOFs from a pose that has ParametersSets.
/// @details For each ParametersSet, for each Parameters, for each perturbable
/// RealValuedParameter (those with can_be_perturbed()==true), creates a ParametricDOFInfo.
/// For parameters marked as global_for_parameters_set, only one DOF is created
/// (associated with the first Parameters object).
void
enumerate_parametric_dofs(
	pose::Pose const & pose,
	utility::vector1< ParametricDOFInfo > & dof_infos
) {
	using namespace core::conformation::parametric;

	dof_infos.clear();

	core::Size const n_param_sets( pose.conformation().n_parameters_sets() );
	if ( n_param_sets == 0 ) {
		TR.Warning << "No ParametersSets found in the pose conformation. No parametric DOFs to enumerate." << std::endl;
		return;
	}

	for ( core::Size ps_index = 1; ps_index <= n_param_sets; ++ps_index ) {

		// Const access to the ParametersSet.
		ParametersSetCOP params_set( pose.conformation().parameters_set( ps_index ) );
		runtime_assert_string_msg( params_set != nullptr,
			"Error in core::optimization::enumerate_parametric_dofs(): Null ParametersSet encountered." );

		core::Size const n_params( params_set->n_parameters() );

		// Keep track of which global parameters we have already added (by parameter index within a Parameters object).
		// We only want to add each global parameter once per ParametersSet.
		std::set< core::Size > global_params_added;

		for ( core::Size p_index = 1; p_index <= n_params; ++p_index ) {

			ParametersCOP params( params_set->parameters( p_index ) );
			runtime_assert_string_msg( params != nullptr,
				"Error in core::optimization::enumerate_parametric_dofs(): Null Parameters encountered." );

			core::Size const n_individual_params( params->num_parameters() );

			for ( core::Size ip = 1; ip <= n_individual_params; ++ip ) {

				ParameterCOP param( params->parameter_cop( ip ) );
				if ( param == nullptr ) continue;

				// Only consider RealValuedParameters.
				RealValuedParameterCOP real_param(
					utility::pointer::dynamic_pointer_cast< RealValuedParameter const >( param ) );
				if ( real_param == nullptr ) continue;

				// Only consider parameters that can be perturbed.
				if ( !real_param->can_be_perturbed() ) continue;

				// For global parameters, only add once (associated with the first Parameters object).
				if ( real_param->global_for_parameters_set() ) {
					if ( global_params_added.count( ip ) > 0 ) continue;
					global_params_added.insert( ip );
				}

				ParametricDOFInfo info;
				info.params_set_index = ps_index;
				info.params_index = p_index;
				info.param_enum = ip;
				info.helix_start_resid = params->first_residue_index();
				info.helix_end_resid = params->last_residue_index();
				info.is_global = real_param->global_for_parameters_set();
				info.scale_factor = 1.0;
				info.param_name = real_param->parameter_name();

				dof_infos.push_back( info );

				TR.Debug << "Enumerated parametric DOF: params_set=" << ps_index
					<< " params=" << p_index
					<< " param=" << ip
					<< " name=" << real_param->parameter_name()
					<< " global=" << ( info.is_global ? "true" : "false" )
					<< " residues=" << info.helix_start_resid << "-" << info.helix_end_resid
					<< std::endl;
			}
		}
	}

	TR << "Enumerated " << dof_infos.size() << " parametric DOFs from " << n_param_sets << " ParametersSet(s)." << std::endl;
}

/// @brief Get the set of residue indices that are under parametric control.
/// @details Reads the residue lists from all Parameters objects in all ParametersSets.
std::set< Size >
get_parametric_residues(
	pose::Pose const & pose
) {
	using namespace core::conformation::parametric;

	std::set< Size > parametric_resids;

	core::Size const n_param_sets( pose.conformation().n_parameters_sets() );
	for ( core::Size ps_index = 1; ps_index <= n_param_sets; ++ps_index ) {

		ParametersSetCOP params_set( pose.conformation().parameters_set( ps_index ) );
		if ( params_set == nullptr ) continue;

		core::Size const n_params( params_set->n_parameters() );
		for ( core::Size p_index = 1; p_index <= n_params; ++p_index ) {

			ParametersCOP params( params_set->parameters( p_index ) );
			if ( params == nullptr ) continue;

			core::Size const first_res( params->first_residue_index() );
			core::Size const last_res( params->last_residue_index() );
			for ( core::Size res = first_res; res <= last_res; ++res ) {
				parametric_resids.insert( res );
			}
		}
	}

	return parametric_resids;
}

/// @brief Get the current value of a parametric DOF from the pose.
Real
get_parametric_dof_value(
	pose::Pose const & pose,
	ParametricDOFInfo const & info
) {
	using namespace core::conformation::parametric;

	runtime_assert_string_msg( info.params_set_index >= 1 && info.params_set_index <= pose.conformation().n_parameters_sets(),
		"Error in core::optimization::get_parametric_dof_value(): params_set_index out of range." );

	ParametersSetCOP params_set( pose.conformation().parameters_set( info.params_set_index ) );
	runtime_assert_string_msg( params_set != nullptr,
		"Error in core::optimization::get_parametric_dof_value(): Null ParametersSet." );

	runtime_assert_string_msg( info.params_index >= 1 && info.params_index <= params_set->n_parameters(),
		"Error in core::optimization::get_parametric_dof_value(): params_index out of range." );

	ParametersCOP params( params_set->parameters( info.params_index ) );
	runtime_assert_string_msg( params != nullptr,
		"Error in core::optimization::get_parametric_dof_value(): Null Parameters." );

	runtime_assert_string_msg( info.param_enum >= 1 && info.param_enum <= params->num_parameters(),
		"Error in core::optimization::get_parametric_dof_value(): param_enum out of range." );

	ParameterCOP param( params->parameter_cop( info.param_enum ) );
	RealValuedParameterCOP real_param(
		utility::pointer::dynamic_pointer_cast< RealValuedParameter const >( param ) );
	runtime_assert_string_msg( real_param != nullptr,
		"Error in core::optimization::get_parametric_dof_value(): Parameter is not a RealValuedParameter." );

	return real_param->value();
}

/// @brief Set a parametric DOF value in the pose's ParametersSet.
void
set_parametric_dof_value(
	pose::Pose & pose,
	ParametricDOFInfo const & info,
	Real value
) {
	using namespace core::conformation::parametric;

	runtime_assert_string_msg( info.params_set_index >= 1 && info.params_set_index <= pose.conformation().n_parameters_sets(),
		"Error in core::optimization::set_parametric_dof_value(): params_set_index out of range." );

	ParametersSetOP params_set( pose.conformation().parameters_set( info.params_set_index ) );
	runtime_assert_string_msg( params_set != nullptr,
		"Error in core::optimization::set_parametric_dof_value(): Null ParametersSet." );

	runtime_assert_string_msg( info.params_index >= 1 && info.params_index <= params_set->n_parameters(),
		"Error in core::optimization::set_parametric_dof_value(): params_index out of range." );

	ParametersOP params( params_set->parameters( info.params_index ) );
	runtime_assert_string_msg( params != nullptr,
		"Error in core::optimization::set_parametric_dof_value(): Null Parameters." );

	runtime_assert_string_msg( info.param_enum >= 1 && info.param_enum <= params->num_parameters(),
		"Error in core::optimization::set_parametric_dof_value(): param_enum out of range." );

	ParameterOP param( params->parameter_op( info.param_enum ) );
	RealValuedParameterOP real_param(
		utility::pointer::dynamic_pointer_cast< RealValuedParameter >( param ) );
	runtime_assert_string_msg( real_param != nullptr,
		"Error in core::optimization::set_parametric_dof_value(): Parameter is not a RealValuedParameter." );

	real_param->set_value( value, true /*ignore_use_degrees*/ );
}

} // namespace optimization
} // namespace core
