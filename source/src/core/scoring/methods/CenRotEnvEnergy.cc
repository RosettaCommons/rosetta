// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CenRotEnvEnergy.cc
/// @brief  CenRot version of centroid env
/// @author Yuan Liu


// Unit headers
#include <core/scoring/methods/CenRotEnvEnergy.hh>
#include <core/scoring/methods/CenRotEnvEnergyCreator.hh>

// Package headers
#include <core/scoring/CenRotEnvPairPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the CenRotEnvEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
CenRotEnvEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new CenRotEnvEnergy;
}

ScoreTypes
CenRotEnvEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( cen_rot_env );
	return sts;
}


/// c-tor
CenRotEnvEnergy::CenRotEnvEnergy() :
	parent( new CenRotEnvEnergyCreator ),
	potential_( ScoringManager::get_instance()->get_CenRotEnvPairPotential() ) {}


/// clone
EnergyMethodOP
CenRotEnvEnergy::clone() const {
	return new CenRotEnvEnergy;
}


/// setup
void
CenRotEnvEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const {
	// compute interpolated number of neighbors at various distance cutoffs
	pose.update_residue_neighbors();
	potential_.compute_centroid_environment( pose );
	potential_.compute_dcentroid_environment( pose );
}

void
CenRotEnvEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const {
	// compute interpolated number of neighbors at various distance cutoffs
	pose.update_residue_neighbors();
	potential_.compute_centroid_environment( pose );
	potential_.compute_dcentroid_environment( pose );
}

/// ENV SCORE
void
CenRotEnvEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const {
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) return;
	if ( rsd.aa() > core::chemical::num_canonical_aas ) return;

	/// accumulate total energies
	Real env_score( 0.0 );
	potential_.evaluate_cen_rot_env_score( pose, rsd, env_score );

	emap[ cen_rot_env ] += env_score;
} // residue_energy


void
CenRotEnvEnergy::eval_residue_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const {
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) return;
	if ( rsd.aa() > core::chemical::num_canonical_aas ) return;

	Real weight_env = weights[ cen_rot_env ];

	numeric::xyzVector<Real> f2_cen, f2_cb;
	potential_.evaluate_cen_rot_env_deriv( pose, rsd, f2_cen, f2_cb);

	Vector const cen_x = rsd.atom( rsd.nbr_atom()+1 ).xyz();
	Vector const cen_y = -f2_cen + cen_x;
	Vector const f1_cen( cen_x.cross( cen_y ) );

	Vector const cb_x = rsd.atom( rsd.nbr_atom() ).xyz();
	Vector const cb_y = -f2_cb + cb_x;
	Vector const f1_cb( cb_x.cross( cb_y ) );

	atom_derivs[ rsd.nbr_atom()+1 ].f1() += weight_env * f1_cen;
	atom_derivs[ rsd.nbr_atom()+1 ].f2() += weight_env * f2_cen;
	atom_derivs[ rsd.nbr_atom()].f1() += weight_env * f1_cb;
	atom_derivs[ rsd.nbr_atom()].f2() += weight_env * f2_cb;
}


void
CenRotEnvEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap &
) const {
	potential_.finalize( pose );
}

core::Size
CenRotEnvEnergy::version() const {
	return 1; // Initial versioning
}

}
}
}
