// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CenRotPairEnergy.cc
/// @brief  CenRot version of cen pair
/// @author Yuan Liu


// Unit headers
#include <core/scoring/methods/CenRotPairEnergy.hh>
#include <core/scoring/methods/CenRotPairEnergyCreator.hh>

// Package headers
#include <core/scoring/CenRotEnvPairPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/DerivVectorPair.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>

#include <numeric/xyzVector.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/angle_deriv.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace methods {

methods::EnergyMethodOP
CenRotPairEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new methods::CenRotPairEnergy );
}

/// @brief Return the set of score types claimed by the EnergyMethod
/// this EnergyMethodCreator creates in its create_energy_method() function
ScoreTypes
CenRotPairEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( cen_rot_pair );
	sts.push_back( cen_rot_pair_ang );
	sts.push_back( cen_rot_pair_dih );
	return sts;
}

/// c-tor
CenRotPairEnergy::CenRotPairEnergy() :
	parent( methods::EnergyMethodCreatorOP( new CenRotPairEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_CenRotEnvPairPotential() )
{}


/// clone
EnergyMethodOP
CenRotPairEnergy::clone() const {
	return EnergyMethodOP( new CenRotPairEnergy() );
}


/////////////////////////////////////////////////////////////////////////////
// score
/////////////////////////////////////////////////////////////////////////////
void
CenRotPairEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	// compute interpolated number of neighbors at various distance cutoffs
	pose.update_residue_neighbors();
}

void
CenRotPairEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd1.has_variant_type( core::chemical::REPLONLY ) || rsd2.has_variant_type( core::chemical::REPLONLY ) ) return;
	if ( rsd1.aa() > core::chemical::num_canonical_aas || rsd2.aa() > core::chemical::num_canonical_aas ) return;

	conformation::Atom const & cen1 ( rsd1.atom( rsd1.nbr_atom()+1 ));
	conformation::Atom const & cen2 ( rsd2.atom( rsd2.nbr_atom()+1 ));
	Real const cendist = cen1.xyz().distance( cen2.xyz() );

	// important: the real cutoff is: atomic_interaction_cutoff + nbr_rad1 + nbr_rad2 !!!
	// maybe not ... use bigger cutoff to be safe 
	// check the boundary directly
	if (cendist<1e-4 || cendist>12.0) return;

	// accumulate total energies
	Real pair_score( 0.0 );
	potential_.evaluate_cen_rot_pair_score( rsd1, rsd2, cendist, pair_score );
	Real pair_ang1_score( 0.0 );
	Real pair_ang2_score( 0.0 );
	Real pair_dih_score( 0.0 );
	potential_.evaluate_cen_rot_pair_orientation_score( rsd1, rsd2, cendist, 
		pair_ang1_score, pair_ang2_score, pair_dih_score );

	emap[ cen_rot_pair ] += pair_score;
	emap[ cen_rot_pair_ang ] += pair_ang1_score + pair_ang2_score;
	emap[ cen_rot_pair_dih ] += pair_dih_score;
}

void
CenRotPairEnergy::eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const &,
		pose::Pose const &,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const {
	if ( rsd1.has_variant_type( core::chemical::REPLONLY ) || rsd2.has_variant_type( core::chemical::REPLONLY ) ) return;
	if ( rsd1.aa() > core::chemical::num_canonical_aas || rsd2.aa() > core::chemical::num_canonical_aas ) return;

	Real weight = weights[ cen_rot_pair ];

	/// assumes centroids are being used, and it is next to nbr atom
	Vector const atom_x = rsd1.atom(rsd1.nbr_atom()+1).xyz();
	Vector const atom_y = rsd2.atom(rsd2.nbr_atom()+1).xyz();
	Real const cendist = atom_x.distance( atom_y );

	if (cendist<1e-4 || cendist>12.0) return;

	Real dpair_dr(0.0);
	potential_.evaluate_cen_rot_pair_deriv( rsd1, rsd2, cendist, dpair_dr );

	Vector f1( atom_x.cross( atom_y ) );
	Vector f2( atom_x - atom_y );
	dpair_dr /= cendist;

	f1 *= weight * dpair_dr;
	f2 *= weight * dpair_dr;

	r1_atom_derivs[ rsd1.nbr_atom()+1 ].f1() += f1;
	r1_atom_derivs[ rsd1.nbr_atom()+1 ].f2() += f2;
	r2_atom_derivs[ rsd2.nbr_atom()+1 ].f1() -= f1;
	r2_atom_derivs[ rsd2.nbr_atom()+1 ].f2() -= f2;

	//orient
	if (rsd1.aa()==chemical::aa_gly || rsd1.aa()==chemical::aa_ala || rsd2.aa()==chemical::aa_gly || rsd2.aa()==chemical::aa_ala) return;

	using namespace numeric::deriv;

	weight = weights[ cen_rot_pair_ang ];
	Real dummy_dih(0.0); //no dihedral yet
	Real dE_dr(0.0);

	//angle 1
	Real theta1(0.0), dE_dtheta1;
	angle_p1_deriv(rsd1.atom(rsd1.nbr_atom()).xyz(), atom_x, atom_y, theta1, f1, f2);
	potential_.evaluate_cen_rot_pair_orientation_deriv(rsd1, rsd2, cendist, theta1, dE_dr, dE_dtheta1, dummy_dih);
	dE_dtheta1 *= weight;
	r1_atom_derivs[rsd1.nbr_atom()].f1() += f1 * dE_dtheta1;
	r1_atom_derivs[rsd1.nbr_atom()].f2() += f2 * dE_dtheta1;
	angle_p2_deriv(rsd1.atom(rsd1.nbr_atom()).xyz(), atom_x, atom_y, theta1, f1, f2);
	r1_atom_derivs[rsd1.nbr_atom()+1].f1() += f1 * dE_dtheta1;
	r1_atom_derivs[rsd1.nbr_atom()+1].f2() += f2 * dE_dtheta1;
	angle_p1_deriv(atom_y, atom_x, rsd1.atom(rsd1.nbr_atom()).xyz(), theta1, f1, f2);
	r2_atom_derivs[rsd2.nbr_atom()+1].f1() += f1 * dE_dtheta1;
	r2_atom_derivs[rsd2.nbr_atom()+1].f2() += f2 * dE_dtheta1;
	Real dis;
	distance_f1_f2_deriv(atom_x, atom_y, dis, f1, f2);
debug_assert( dis==cendist );
	dE_dr *= weight;
	r1_atom_derivs[rsd1.nbr_atom()+1].f1() += f1 * dE_dr;
	r1_atom_derivs[rsd1.nbr_atom()+1].f2() += f2 * dE_dr;
	r2_atom_derivs[rsd2.nbr_atom()+1].f1() -= f1 * dE_dr;
	r2_atom_derivs[rsd2.nbr_atom()+1].f2() -= f2 * dE_dr;

	//angle 2
	Real theta2(0.0), dE_dtheta2;
	angle_p1_deriv(rsd2.atom(rsd2.nbr_atom()).xyz(), atom_y, atom_x, theta2, f1, f2);
	potential_.evaluate_cen_rot_pair_orientation_deriv(rsd2, rsd1, cendist, theta2, dE_dr, dE_dtheta2, dummy_dih);
	dE_dtheta2 *= weight;
	r2_atom_derivs[rsd2.nbr_atom()].f1() += f1 * dE_dtheta2;
	r2_atom_derivs[rsd2.nbr_atom()].f2() += f2 * dE_dtheta2;
	angle_p1_deriv(atom_x, atom_y, rsd2.atom(rsd2.nbr_atom()).xyz(), theta2, f1, f2);
	r1_atom_derivs[rsd1.nbr_atom()+1].f1() += f1 * dE_dtheta2;
	r1_atom_derivs[rsd1.nbr_atom()+1].f2() += f2 * dE_dtheta2;
	angle_p2_deriv(rsd2.atom(rsd2.nbr_atom()).xyz(), atom_y, atom_x, theta2, f1, f2);
	r2_atom_derivs[rsd2.nbr_atom()+1].f1() += f1 * dE_dtheta2;
	r2_atom_derivs[rsd2.nbr_atom()+1].f2() += f2 * dE_dtheta2;

	distance_f1_f2_deriv(atom_x, atom_y, dis, f1, f2);
debug_assert( dis==cendist );
	dE_dr *= weight;
	r1_atom_derivs[rsd1.nbr_atom()+1].f1() += f1 * dE_dr;
	r1_atom_derivs[rsd1.nbr_atom()+1].f2() += f2 * dE_dr;
	r2_atom_derivs[rsd2.nbr_atom()+1].f1() -= f1 * dE_dr;
	r2_atom_derivs[rsd2.nbr_atom()+1].f2() -= f2 * dE_dr;
}


void
CenRotPairEnergy::finalize_total_energy(
	pose::Pose &,
	ScoreFunction const &,
	EnergyMap &
) const {
}

/// @brief CenRotPairEnergy distance cutoff
Distance
CenRotPairEnergy::atomic_interaction_cutoff() const {
	return 12;
	//return 6.0;
}

core::Size
CenRotPairEnergy::version() const {
	return 1; // Initial versioning
}

}
}
}
