// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/energy_methods/FaMPEnvEnergy.hh
///
/// @brief  LK-Type Membrane Environment Energy
///
/// @author  Patrick Barth (Original)
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_FaMPEnvEnergy_cc
#define INCLUDED_core_scoring_membrane_FaMPEnvEnergy_cc

// Unit headers
#include <core/energy_methods/FaMPEnvEnergy.hh>
#include <core/energy_methods/FaMPEnvEnergyCreator.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/MembraneGeometry.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/memb_etable/MembEtable.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

// Utility Headers
#include <ObjexxFCL/FArray1.fwd.hh>
#include <utility>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers

static basic::Tracer TR( "core.energy_methods.FaMPEnvEnergy" );

namespace core {
namespace energy_methods {

using namespace core::scoring::methods;

/// @brief Return a fresh instance of the energy method
core::scoring::methods::EnergyMethodOP
FaMPEnvEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< FaMPEnvEnergy >( ( core::scoring::ScoringManager::get_instance()->memb_etable( options.etable_type() ) ) );
}

/// @brief Return Score Types Required for Method
core::scoring::ScoreTypes
FaMPEnvEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( FaMPEnv );
	return sts;
}

/// @brief Construct Energy Method from Etable
FaMPEnvEnergy::FaMPEnvEnergy( core::scoring::etable::MembEtableCAP memb_etable_in ) :
	parent( utility::pointer::make_shared< FaMPEnvEnergyCreator >() ),
	memb_etable_(std::move( memb_etable_in )),
	lk_dgrefce_( memb_etable_.lock()->lk_dgrefce() ), // FIXME: check lock() success
	memb_lk_dgrefce_( memb_etable_.lock()->memb_lk_dgrefce() ) // FIXME: check lock() success
{}

/// @brief Clone Energy Method
core::scoring::methods::EnergyMethodOP
FaMPEnvEnergy::clone() const {
	return utility::pointer::make_shared< FaMPEnvEnergy >( *this );
}

/// @brief Compute Per-Residue Energies
void
FaMPEnvEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	core::scoring::EnergyMap & emap
) const {

	for ( Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end; ++i ) {
		emap[ core::scoring::FaMPEnv ] += eval_fa_mbenv( rsd.atom(i), fa_proj_[ rsd.seqpos() ][ i ] );
	}
}

/// @brief Fianlzie Total Per-Residue Energies
void
FaMPEnvEnergy::finalize_total_energy(
	pose::Pose &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const {

	emap[ core::scoring::FaMPEnv ] += 0.0;
}

/// @brief Setup for Computing Derivatives
void
FaMPEnvEnergy::setup_for_derivatives(
	pose::Pose & pose,
	core::scoring::ScoreFunction const & scfxn
) const {

	using namespace core::scoring::membrane;

	// Setup Projection and Derivatives
	init( pose );

	// Set membrane env weight from map
	fa_mbenv_weight_ = scfxn.weights()[ core::scoring::FaMPEnv ];
}

/// @brief Evaluate Per-Atom Derivatives
void
FaMPEnvEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const &,
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


	Real cp_weight = 1.0;

	// Initialize f1, f2
	Vector f1( 0.0 ), f2( 0.0 );


	Real const deriv = memb_lk_dgrefce_( rsd1.atom( m ).type() ) - lk_dgrefce_( rsd1.atom( m ).type() );

	//F1 and F2 are defined in H. Abe, W. Braun, T. Noguti, N. Go, Computers & Chemistry 1984

	f1 = fa_f1_[ rsd1.seqpos() ][ m ];
	f2 = fa_f2_[ rsd1.seqpos() ][ m ];


	F1 += fa_mbenv_weight_* deriv * cp_weight * f1;
	F2 += fa_mbenv_weight_ * deriv * cp_weight * f2;





}

/// @brief Fa_MbenvEnergy is context independent
void
FaMPEnvEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const {}

/// @brief Setup Method for initial scoring
void
FaMPEnvEnergy::setup_for_scoring(
	pose::Pose & pose,
	core::scoring::ScoreFunction const &
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
	core::conformation::membrane::MembraneGeometryCOP mp_geometry( conf.membrane_info()->membrane_geometry() );


	// For convenience - grab nres
	Real nres = pose.size();

	for ( Size i = 1; i <= nres; ++i ) {
		for ( Size j = 1, j_end = pose.residue( i ).nheavyatoms(); j <= j_end; ++j ) {

			Vector const xyz( pose.residue( i ).xyz( j ) );

			//Compute transition function
			fa_proj_[i][j] = mp_geometry->f_transition( conf, i, j );


			// Compute derivatives
			fa_f1_[i][j] = mp_geometry->f_transition_f1( conf, i, j );
			fa_f2_[i][j] = mp_geometry->f_transition_f2( conf, i, j );

		}
	}
}


/// @brief Setup Data Members for Fullatom Info
void
FaMPEnvEnergy::setup_for_fullatom( pose::Pose & pose ) const {

	core::Real nres = pose.size();

	fa_proj_.resize( (core::Size)nres );
	fa_f1_.resize( (core::Size)nres );
	fa_f2_.resize( (core::Size)nres );

	static Size const MAX_AMINOACID_SIZE = 15;

	for ( Size i = 1; i <= nres; ++i ) {

		Size const max_size = std::max( MAX_AMINOACID_SIZE, pose.residue( i ).nheavyatoms() );

		fa_proj_[i].resize( max_size );
		fa_f1_[i].resize( max_size );
		fa_f2_[i].resize( max_size );

		for ( Size j = 1; j <= max_size; ++j ) {

			fa_proj_[i][j] = 0.0;
			fa_f1_[i][j].assign(0.0,0.0,0.0);
			fa_f2_[i][j].assign(0.0,0.0,0.0);

		}
	}
}


} // scoring
} // core

#endif // INCLUDED_core_energy_methods_FaMPEnvEnergy_cc
