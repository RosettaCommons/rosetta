// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/FaMPEnvEnergy.hh
///
/// @brief  LK-Type Membrane Environment Energy
/// @details Last Modified: 5/13/14
///
/// @author  Patrick Barth (Original)
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_FaMPEnvEnergy_cc
#define INCLUDED_core_scoring_membrane_FaMPEnvEnergy_cc

// Unit headers
#include <core/scoring/membrane/FaMPEnvEnergy.hh>
#include <core/scoring/membrane/FaMPEnvEnergyCreator.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>

// Project Headers
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/memb_etable/MembEtable.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

// Utility Headers
#include <ObjexxFCL/FArray1.fwd.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.membrane.FaMPEnvEnergy" );

namespace core {
namespace scoring {
namespace membrane {

using namespace core::scoring::methods;

/// @brief Return a fresh instance of the energy method
methods::EnergyMethodOP
FaMPEnvEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new FaMPEnvEnergy( ( ScoringManager::get_instance()->memb_etable( options.etable_type() ) ) ) );
}

/// @brief Return Score Types Required for Method
ScoreTypes
FaMPEnvEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( FaMPEnv );
	return sts;
}

/// @brief Construct Energy Method from Etable
FaMPEnvEnergy::FaMPEnvEnergy( etable::MembEtableCAP memb_etable_in ) :
	parent( EnergyMethodCreatorOP( new FaMPEnvEnergyCreator ) ),
	memb_etable_( memb_etable_in ),
	lk_dgrefce_( memb_etable_.lock()->lk_dgrefce() ), // FIXME: check lock() success
	memb_lk_dgrefce_( memb_etable_.lock()->memb_lk_dgrefce() ) // FIXME: check lock() success
{}

/// @brief Clone Energy Method
EnergyMethodOP
FaMPEnvEnergy::clone() const {
	return EnergyMethodOP( new FaMPEnvEnergy( *this ) );
}

/// @brief Compute Per-Residue Energies
void
FaMPEnvEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const {

	for ( Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end; ++i ) {
		emap[ FaMPEnv ] += eval_fa_mbenv( rsd.atom(i), fa_proj_[ rsd.seqpos() ][ i ] );
	}
}

/// @brief Fianlzie Total Per-Residue Energies
void
FaMPEnvEnergy::finalize_total_energy(
	pose::Pose &,
	ScoreFunction const &,
	EnergyMap & emap
) const {

	emap[ FaMPEnv ] += 0.0;
}

/// @brief Setup for Computing Derivatives
void
FaMPEnvEnergy::setup_for_derivatives(
	pose::Pose & pose,
	ScoreFunction const & scfxn
) const {

	using namespace core::scoring::membrane;

	// Setup Projection and Derivatives
	init( pose );

	// Set membrane env weight from map
	fa_mbenv_weight_ = scfxn.weights()[ FaMPEnv ];
}

/// @brief Evaluate Per-Atom Derivatives
void
FaMPEnvEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const &,
	Vector & F1,
	Vector & F2
) const {

	// Grab residue and atom numbers
	Size const i( atom_id.rsd() );
	Size const m( atom_id.atomno() );

	// Get the actual residue and atom
	conformation::Residue const & rsd1( pose.residue( i ) );
	if ( m > rsd1.nheavyatoms() ) return;
	Vector const heavy_atom_i( rsd1.xyz( m ) );

	// Get Membrane Center from Pose
	core::conformation::Conformation const & conf( pose.conformation() );
	core::Vector center = conf.membrane_info()->membrane_center(conf);
	Real cp_weight = 1.0;

	// Initialize f1, f2
	Vector f1( 0.0 ), f2( 0.0 );

	// Grab derivative of proj
	core::Real fa_proj_deriv = fa_proj_deriv_[ rsd1.seqpos() ][ m ];
	core::Vector fa_proj_coord = fa_proj_coord_[ rsd1.seqpos() ][ m ];

	Real const deriv = memb_lk_dgrefce_( rsd1.atom( m ).type() ) - lk_dgrefce_( rsd1.atom( m ).type() );
	Real dE_dZ_over_r = fa_mbenv_weight_ * deriv * fa_proj_deriv;

	Vector const d_ij = fa_proj_coord - heavy_atom_i;
	Real const d_ij_norm = d_ij.length();
	if ( d_ij_norm == Real(0.0) ) return;

	Real const invd = 1.0 / d_ij_norm;
	f2 = d_ij * invd;
	f1 = fa_proj_coord.cross(heavy_atom_i);
	f1 *= invd;

	if ( dE_dZ_over_r != 0.0 ) {
		F1 += dE_dZ_over_r * cp_weight * f1;
		F2 += dE_dZ_over_r * cp_weight * f2;
	}
}

/// @brief Fa_MbenvEnergy is context independent
void
FaMPEnvEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const {}

/// @brief Setup Method for initial scoring
void
FaMPEnvEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const &
) const {

	// Setup Projection and Derivatives
	init( pose );
}

/// @brief Evaluate Per-Atom Env term
Real
FaMPEnvEnergy::eval_fa_mbenv(
	conformation::Atom const & atom1,
	Real const & f1
) const {

	Real score( 0.0 );
	score = ( 1 - f1 ) * ( memb_lk_dgrefce_( atom1.type() ) - lk_dgrefce_( atom1.type() ) );
	return score;
}

/// @brief Versioning
core::Size
FaMPEnvEnergy::version() const {
	return 2;
}

void
FaMPEnvEnergy::init( pose::Pose & pose ) const {

	// Alloc appropriate Farrays
	setup_for_fullatom( pose );

	// Grab membrane center/normal from pose conf
	core::conformation::Conformation const & conf( pose.conformation() );
	core::Vector center = conf.membrane_info()->membrane_center(conf);
	core::Vector normal = conf.membrane_info()->membrane_normal(conf);

	core::Real thickness = conf.membrane_info()->membrane_thickness();
	core::Real steepness = conf.membrane_info()->membrane_steepness();

	// For convenience - grab nres
	Real nres = pose.total_residue();

	for ( Size i = 1; i <= nres; ++i ) {
		for ( Size j = 1, j_end = pose.residue( i ).nheavyatoms(); j <= j_end; ++j ) {

			Vector const xyz( pose.residue( i ).xyz( j ) );

			// Compute Standard Z Position
			fa_z_position_[i][j] = conf.membrane_info()->atom_z_position( conf, i, j );

			// Compute Fa Projection
			fa_proj_[i][j] = compute_fa_proj( fa_z_position_[i][j], thickness, steepness );

			// Compute derivatives
			fa_proj_deriv_[i][j] = compute_fa_deriv( fa_z_position_[i][j], thickness, steepness );

			// Compute Projected coordinate
			fa_proj_coord_[i][j] = compute_fa_proj_coord( fa_z_position_[i][j], xyz, center, normal );

		}
	}
}

/// @brief Helper Method - Compute Fa Proj
core::Real
FaMPEnvEnergy::compute_fa_proj(
	core::Real z_position,
	core::Real thickness,
	core::Real steepness
) const {

	Real internal_product(0), z(0), zn(0);
	internal_product = std::abs( z_position );
	z = internal_product;
	z /= thickness;
	zn = std::pow( z, steepness );
	Real result = zn/(1 + zn);

	return result;
}

/// @brief Helper Method - Compute Fa Derivatives
core::Real
FaMPEnvEnergy::compute_fa_deriv(
	core::Real z_position,
	core::Real thickness,
	core::Real steepness
) const {

	Real internal_product(0), z(0), zn(0), znm1(0);
	internal_product = std::abs( z_position );
	z = internal_product;
	z /= thickness;
	zn = std::pow( z, steepness );
	znm1 = std::pow( z, steepness-1 );
	Real deriv = steepness * znm1 * std::pow((1+zn),-2);

	return ( deriv /= thickness );
}

/// @brief Helper Method - Compute Fa Proj coordinate
core::Vector
FaMPEnvEnergy::compute_fa_proj_coord(
	core::Real z_position,
	core::Vector xyz,
	core::Vector center,
	core::Vector normal
) const {
	Vector proj_i = center + z_position * normal;
	Vector i_ip = proj_i - xyz;
	return ( center - i_ip );
}


/// @brief Setup Data Members for Fullatom Info
void
FaMPEnvEnergy::setup_for_fullatom( pose::Pose & pose ) const {

	core::Real nres = pose.total_residue();

	fa_proj_.resize( (core::Size)nres );
	fa_proj_coord_.resize( (core::Size)nres );
	fa_proj_deriv_.resize( (core::Size)nres );
	fa_z_position_.resize( (core::Size)nres );

	static Size const MAX_AMINOACID_SIZE = 15;

	for ( Size i = 1; i <= nres; ++i ) {

		Size const max_size = std::max( MAX_AMINOACID_SIZE, pose.residue( i ).nheavyatoms() );

		fa_proj_[i].resize( max_size );
		fa_proj_coord_[i].resize( max_size );
		fa_proj_deriv_[i].resize( max_size );
		fa_z_position_[i].resize( max_size );

		for ( Size j = 1; j <= max_size; ++j ) {

			fa_proj_[i][j] = 0.0;
			fa_proj_coord_[i][j].assign(0.0,0.0,0.0);
			fa_proj_deriv_[i][j] = 0.0;
			fa_z_position_[i][j] = 0.0;

		}
	}
}


} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_methods_FaMPEnvEnergy_cc
