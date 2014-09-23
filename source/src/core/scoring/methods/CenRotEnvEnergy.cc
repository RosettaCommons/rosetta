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
	sts.push_back( cen_rot_cbeta );
	return sts;
}


/// c-tor
CenRotEnvEnergy::CenRotEnvEnergy() :
	parent( methods::EnergyMethodCreatorOP( new CenRotEnvEnergyCreator ) ),
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
	Real cbeta6_score( 0.0 );
	Real cbeta12_score( 0.0 );
	potential_.evaluate_cen_rot_env_and_cbeta_score( pose, rsd, env_score, cbeta6_score, cbeta12_score );

	emap[ cen_rot_env ] += env_score;
	emap[ cen_rot_cbeta ] += cbeta6_score + cbeta12_score;
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
	Real weight_cbeta = weights[ cen_rot_cbeta ];

	numeric::xyzVector<Real> f2_cen_env, f2_cen_cb6, f2_cen_cb12;
	numeric::xyzVector<Real> f2_cb_env, f2_cb_cb6, f2_cb_cb12;
	potential_.evaluate_cen_rot_env_and_cbeta_deriv( pose, rsd,
		f2_cen_env, f2_cen_cb6, f2_cen_cb12,
		f2_cb_env, f2_cb_cb6, f2_cb_cb12);

	//force on CEN
	Vector const cen_x = rsd.atom("CEN").xyz();
	Vector cen_y = -f2_cen_env + cen_x;
	Vector const f1_cen_env( cen_x.cross( cen_y ) );

	Vector const f2_cen_cbeta( f2_cen_cb6 + f2_cen_cb12 );
	cen_y = -f2_cen_cbeta + cen_x;
	Vector const f1_cen_cbeta( cen_x.cross( cen_y ) );

	//force on CB
	Vector const cb_x = rsd.atom( rsd.nbr_atom() ).xyz();
	Vector cb_y = -f2_cb_env + cb_x;
	Vector const f1_cb_env( cb_x.cross( cb_y ) );
	
	Vector const f2_cb_cbeta( f2_cb_cb6 + f2_cb_cb12 );
	cb_y = -f2_cb_cbeta + cb_x;
	Vector const f1_cb_cbeta( cb_x.cross( cb_y ) );

	atom_derivs[ rsd.atom_index("CEN") ].f1() += weight_env * f1_cen_env + weight_cbeta * f1_cen_cbeta;
	atom_derivs[ rsd.atom_index("CEN") ].f2() += weight_env * f2_cen_env + weight_cbeta * f2_cen_cbeta;
	atom_derivs[ rsd.nbr_atom()].f1() += weight_env * f1_cb_env + weight_cbeta * f1_cb_cbeta;
	atom_derivs[ rsd.nbr_atom()].f2() += weight_env * f2_cb_env + weight_cbeta * f2_cb_cbeta;
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
