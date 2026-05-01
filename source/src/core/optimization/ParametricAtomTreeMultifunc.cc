// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/ParametricAtomTreeMultifunc.cc
/// @brief  Multifunc for co-minimizing parametric and atom-tree DOFs.
/// @author Andy Watkins

#include <core/optimization/ParametricAtomTreeMultifunc.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/atom_tree_minimize.hh>
#include <core/optimization/parametric_minimize_util.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/DerivVectorPair.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <core/conformation/parametric/RealVectorValuedParameter.hh>
#include <core/conformation/parametric/SizeValuedParameter.hh>
#include <core/conformation/parametric/SizeVectorValuedParameter.hh>

#include <numeric/crick_equations/BundleParams.hh>
#include <numeric/crick_equations/BundleParams_derivatives.hh>

#include <basic/Tracer.hh>
#include <basic/prof.hh>

namespace core {
namespace optimization {

static basic::Tracer TR( "core.optimization.ParametricAtomTreeMultifunc" );

ParametricAtomTreeMultifunc::ParametricAtomTreeMultifunc(
	pose::Pose & pose_in,
	MinimizerMap & min_map_in,
	scoring::ScoreFunction const & scorefxn_in,
	utility::vector1< ParametricDOFInfo > const & parametric_dofs_in,
	bool const deriv_check_in,
	bool const deriv_check_verbose_in
) :
	pose_( pose_in ),
	min_map_( min_map_in ),
	score_function_( scorefxn_in ),
	parametric_dofs_( parametric_dofs_in ),
	n_standard_dofs_( min_map_in.nangles() ),
	deriv_check_( deriv_check_in ),
	deriv_check_verbose_( deriv_check_verbose_in )
{}

ParametricAtomTreeMultifunc::~ParametricAtomTreeMultifunc() = default;

void
ParametricAtomTreeMultifunc::apply_parametric_dofs( Multivec const & vars ) const {
	for ( Size p = 1; p <= parametric_dofs_.size(); ++p ) {
		Real const val = vars[ n_standard_dofs_ + p ];
		set_parametric_dof_value( pose_, parametric_dofs_[p], val );
	}
	rebuild_parametric_backbone( pose_ );
}

Real
ParametricAtomTreeMultifunc::operator()( Multivec const & vars ) const {
	PROF_START( basic::FUNC );

	apply_parametric_dofs( vars );
	min_map_.copy_dofs_to_pose( pose_, vars );
	Real const score = score_function_( pose_ );

	PROF_STOP( basic::FUNC );
	return score;
}

void
ParametricAtomTreeMultifunc::dfunc( Multivec const & vars, Multivec & dE_dvars ) const {
	PROF_START( basic::DFUNC );

	apply_parametric_dofs( vars );

	// Compute standard DOF derivatives via the normal atom_tree_dfunc pipeline.
	// This fills dE_dvars[1..n_standard_dofs_] and populates
	// min_map_.atom_derivatives_ with per-atom DerivVectorPairs.
	atom_tree_dfunc( pose_, min_map_, score_function_, vars, dE_dvars );

	// Extend dE_dvars to include parametric DOF slots
	Size const ntotal = n_standard_dofs_ + parametric_dofs_.size();
	dE_dvars.resize( ntotal, 0.0 );

	// Compute parametric DOF derivatives via the Crick Jacobian chain rule:
	// dE/dParam = Σ_atoms (dE/dXYZ · dXYZ/dParam)
	for ( Size p = 1; p <= parametric_dofs_.size(); ++p ) {
		ParametricDOFInfo const & info = parametric_dofs_[p];
		Real dE_dparam = 0.0;

		// Get current Crick parameter values for Jacobian computation
		conformation::parametric::ParametersSetCOP params_set =
			pose_.conformation().parameters_set( info.params_set_index );
		conformation::parametric::ParametersCOP params =
			params_set->parameters( info.params_index );

		// Extract the Crick parameters needed for derivative computation.
		// We need: r0, omega0, delta_omega0, omega1, z1, delta_omega1_all, plus
		// per-atom r1, delta_omega1, delta_z1.
		// These are stored as parameters in the Parameters object.
		// The exact enum indices depend on the calculator type (BPC or BBPC).
		// For generality, we look up by name.
		Real r0_val = 0, omega0_val = 0, delta_omega0_val = 0;
		Real omega1_val = 0, z1_val = 0, delta_omega1_all_val = 0;
		utility::vector1< Real > r1_vals, delta_omega1_vals, delta_z1_vals;
		//Size residues_per_repeat = 1;

		// Read parameter values by iterating and checking names
		for ( Size i = 1; i <= params->num_parameters(); ++i ) {
			conformation::parametric::ParameterCOP param = params->parameter_cop( i );
			std::string const & name = param->parameter_name();
			if ( name == "r0" ) {
				r0_val = utility::pointer::static_pointer_cast< conformation::parametric::RealValuedParameter const >( param )->value();
			} else if ( name == "omega0" ) {
				omega0_val = utility::pointer::static_pointer_cast< conformation::parametric::RealValuedParameter const >( param )->value();
			} else if ( name == "delta_omega0" ) {
				delta_omega0_val = utility::pointer::static_pointer_cast< conformation::parametric::RealValuedParameter const >( param )->value();
			} else if ( name == "omega1" ) {
				omega1_val = utility::pointer::static_pointer_cast< conformation::parametric::RealValuedParameter const >( param )->value();
			} else if ( name == "z1" ) {
				z1_val = utility::pointer::static_pointer_cast< conformation::parametric::RealValuedParameter const >( param )->value();
			} else if ( name == "delta_omega1" ) {
				delta_omega1_all_val = utility::pointer::static_pointer_cast< conformation::parametric::RealValuedParameter const >( param )->value();
			} else if ( name == "r1_peratom" ) {
				r1_vals = utility::pointer::static_pointer_cast< conformation::parametric::RealVectorValuedParameter const >( param )->value();
			} else if ( name == "delta_omega1_peratom" ) {
				delta_omega1_vals = utility::pointer::static_pointer_cast< conformation::parametric::RealVectorValuedParameter const >( param )->value();
			} else if ( name == "delta_z1_peratom" ) {
				delta_z1_vals = utility::pointer::static_pointer_cast< conformation::parametric::RealVectorValuedParameter const >( param )->value();
			} /*else if ( name == "residues_per_repeat" ) {
				residues_per_repeat = utility::pointer::static_pointer_cast< conformation::parametric::SizeValuedParameter const >( param )->value();
			}*/
		}

		if ( r1_vals.empty() || delta_omega1_vals.empty() || delta_z1_vals.empty() ) continue;

		// Compute omega1 relative to omega0 (as done in generate_atom_positions)
		Real const omega1_relative = omega1_val - omega0_val;

		// Compute t values and iterate over mainchain atoms
		Size const helix_length = info.helix_end_resid - info.helix_start_resid + 1;
		Real t = -static_cast<Real>( helix_length + 2 ) / 2.0;

		Size atom_counter = 0;
		for ( Size resid = info.helix_start_resid; resid <= info.helix_end_resid; ++resid ) {
			Size const n_mainchain = pose_.residue( resid ).n_mainchain_atoms();
			for ( Size iatom = 1; iatom <= n_mainchain; ++iatom ) {
				++atom_counter;
				Size const atom_index_in_repeat = ((atom_counter - 1) % r1_vals.size()) + 1;

				Real const r1 = r1_vals[ atom_index_in_repeat ];
				Real const dw1 = delta_omega1_vals[ atom_index_in_repeat ] + delta_omega1_all_val;
				Real const dz1 = delta_z1_vals[ atom_index_in_repeat ];

				// Compute Crick derivatives for this atom
				bool deriv_failed = false;
				numeric::crick_equations::CrickDerivatives crick_derivs =
					numeric::crick_equations::compute_crick_derivatives(
						t, r0_val, omega0_val, delta_omega0_val,
						r1, omega1_relative, z1_val,
						dw1, dz1, deriv_failed );

				if ( deriv_failed ) continue;

				// Get dE/dXYZ for this atom from the MinimizerMap's atom_derivatives
				Size const real_atomno = pose_.residue_type( resid ).mainchain_atom( iatom );
				scoring::DerivVectorPair const & dvp = min_map_.atom_derivatives( resid )[ real_atomno ];
				numeric::xyzVector< Real > const & dE_dXYZ = dvp.f2();

				// Select the correct derivative component based on which parametric DOF this is
				numeric::xyzVector< Real > dXYZ_dParam( 0.0, 0.0, 0.0 );
				std::string const & dof_name = info.param_name;
				if ( dof_name == "r0" ) {
					dXYZ_dParam = crick_derivs.dXYZ_dr0;
				} else if ( dof_name == "omega0" ) {
					dXYZ_dParam = crick_derivs.dXYZ_domega0;
				} else if ( dof_name == "delta_omega0" ) {
					dXYZ_dParam = crick_derivs.dXYZ_ddelta_omega0;
				} else if ( dof_name == "delta_omega1" ) {
					dXYZ_dParam = crick_derivs.dXYZ_ddelta_omega1;
				} else if ( dof_name == "delta_t" || dof_name == "delta_z0" ) {
					dXYZ_dParam = crick_derivs.dXYZ_ddelta_t;
				}

				dE_dparam += dE_dXYZ.dot( dXYZ_dParam );
			}
			// Advance t by 1 per residue (matching generate_atom_positions logic)
			t += 1.0;
		}

		dE_dvars[ n_standard_dofs_ + p ] = dE_dparam;
	}

	if ( deriv_check_ ) {
		numerical_derivative_check( min_map_, *this, vars, dE_dvars, nullptr, deriv_check_verbose_ );
	}

	PROF_STOP( basic::DFUNC );
}

void
ParametricAtomTreeMultifunc::dump( Multivec const & /*vars*/, Multivec const & /*vars2*/ ) const {
	// Minimal implementation — extend if debugging is needed
}

} // namespace optimization
} // namespace core
