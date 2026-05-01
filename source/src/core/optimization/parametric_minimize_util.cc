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
#include <core/conformation/parametric/RealVectorValuedParameter.hh>
#include <core/conformation/parametric/SizeValuedParameter.hh>
#include <core/conformation/parametric/SizeVectorValuedParameter.hh>
#include <core/conformation/parametric/BooleanValuedParameter.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

#include <numeric/crick_equations/BundleParams.hh>

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

void
rebuild_parametric_backbone(
	pose::Pose & pose
) {
	using namespace core::conformation::parametric;

	for ( core::Size ps = 1; ps <= pose.conformation().n_parameters_sets(); ++ps ) {
		ParametersSetCOP params_set = pose.conformation().parameters_set( ps );
		if ( params_set == nullptr ) continue;
		for ( core::Size p = 1; p <= params_set->n_parameters(); ++p ) {
			ParametersCOP params = params_set->parameters( p );
			if ( params == nullptr ) continue;
			if ( params->n_residue() == 0 ) continue;

			core::Size const start = params->first_residue_index();
			core::Size const end = params->last_residue_index();

			// Rebuild this parametric element's backbone coordinates.
			// We use the helical_bundle utility functions directly: generate new atom
			// positions from current parameter values, then place them.
			// For now, this is a simplified rebuild that re-applies the Crick equations.
			// A more complete version would use the ParametrizationCalculator's build_helix/build_strand.

			// Extract parameters by name for the Crick equation call
			Real r0_val = 0, omega0_val = 0, delta_omega0_val = 0;
			Real omega1_val = 0, z1_val = 0, delta_omega1_all_val = 0;
			Real delta_t_val = 0, epsilon_val = 1.0;
			bool invert = false;
			utility::vector1< Real > r1_vals, delta_omega1_vals, delta_z1_vals;
			core::Size residues_per_repeat = 1, repeating_unit_offset = 0;
			utility::vector1< core::Size > atoms_per_residue_vals;

			for ( core::Size i = 1; i <= params->num_parameters(); ++i ) {
				ParameterCOP param = params->parameter_cop( i );
				std::string const & name = param->parameter_name();
				if ( name == "r0" ) {
					r0_val = utility::pointer::static_pointer_cast< RealValuedParameter const >( param )->value();
				} else if ( name == "omega0" ) {
					omega0_val = utility::pointer::static_pointer_cast< RealValuedParameter const >( param )->value();
				} else if ( name == "delta_omega0" ) {
					delta_omega0_val = utility::pointer::static_pointer_cast< RealValuedParameter const >( param )->value();
				} else if ( name == "omega1" ) {
					omega1_val = utility::pointer::static_pointer_cast< RealValuedParameter const >( param )->value();
				} else if ( name == "z1" ) {
					z1_val = utility::pointer::static_pointer_cast< RealValuedParameter const >( param )->value();
				} else if ( name == "delta_omega1" ) {
					delta_omega1_all_val = utility::pointer::static_pointer_cast< RealValuedParameter const >( param )->value();
				} else if ( name == "delta_t" || name == "delta_z0" ) {
					delta_t_val = utility::pointer::static_pointer_cast< RealValuedParameter const >( param )->value();
				} else if ( name == "epsilon" ) {
					epsilon_val = utility::pointer::static_pointer_cast< RealValuedParameter const >( param )->value();
				} else if ( name == "invert" ) {
					invert = utility::pointer::static_pointer_cast< BooleanValuedParameter const >( param )->value();
				} else if ( name == "r1_peratom" ) {
					r1_vals = utility::pointer::static_pointer_cast< RealVectorValuedParameter const >( param )->value();
				} else if ( name == "delta_omega1_peratom" ) {
					delta_omega1_vals = utility::pointer::static_pointer_cast< RealVectorValuedParameter const >( param )->value();
				} else if ( name == "delta_z1_peratom" ) {
					delta_z1_vals = utility::pointer::static_pointer_cast< RealVectorValuedParameter const >( param )->value();
				} else if ( name == "residues_per_repeat" ) {
					residues_per_repeat = utility::pointer::static_pointer_cast< SizeValuedParameter const >( param )->value();
				} else if ( name == "atoms_per_residue" ) {
					atoms_per_residue_vals = utility::pointer::static_pointer_cast< SizeVectorValuedParameter const >( param )->value();
				}
			}

			if ( r1_vals.empty() ) continue;

			// Compute omega1 relative to omega0 (matching generate_atom_positions convention)
			Real const omega1_relative = omega1_val - omega0_val;

			// Compute atom positions directly via Crick equations (numeric layer, no protocols dependency).
			core::Size const helix_length = end - start + 1;
			Real t = -static_cast<Real>( helix_length + 2 ) / 2.0 + delta_t_val;

			bool rebuild_failed = false;
			core::Size atom_counter = 0;
			for ( core::Size resid = start; resid <= end && !rebuild_failed; ++resid ) {
				core::Size const n_mc = pose.residue( resid ).n_mainchain_atoms();
				for ( core::Size iatom = 1; iatom <= n_mc && !rebuild_failed; ++iatom ) {
					++atom_counter;
					core::Size const idx = ((atom_counter - 1) % r1_vals.size()) + 1;

					Real const r1 = r1_vals[ idx ];
					Real const dw1 = delta_omega1_vals[ idx ] + delta_omega1_all_val;
					Real const dz1 = delta_z1_vals[ idx ];

					bool failed = false;
					numeric::xyzVector< Real > xyz = numeric::crick_equations::XYZ_BUNDLE(
						t, r0_val, omega0_val, delta_omega0_val,
						r1, omega1_relative, z1_val, dw1, dz1, epsilon_val, failed );

					if ( failed ) {
						TR.Warning << "rebuild_parametric_backbone: Crick equation failed for element "
							<< p << " residue " << resid << " atom " << iatom << std::endl;
						rebuild_failed = true;
						break;
					}

					if ( invert ) {
						xyz.x( -xyz.x() );
						xyz.z( -xyz.z() );
					}

					core::Size const real_atomno = pose.residue_type( resid ).mainchain_atom( iatom );
					pose.set_xyz( id::AtomID( real_atomno, resid ), xyz );
				}
				t += 1.0;
			}
		}
	}
}

} // namespace optimization
} // namespace core
