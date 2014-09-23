// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/YHHPlanarityEnergy.hh
/// @brief  Term for chi3 on tyrosine residues to prefer the hydrogen lie in the plane of the ring
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit headers
#include <core/scoring/methods/YHHPlanarityEnergy.hh>
#include <core/scoring/methods/YHHPlanarityEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/constants.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the P_AA_pp_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
YHHPlanarityEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new YHHPlanarityEnergy;
}

ScoreTypes
YHHPlanarityEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( yhh_planarity );
	return sts;
}



/// ctor
YHHPlanarityEnergy::YHHPlanarityEnergy() :
	parent( methods::EnergyMethodCreatorOP( new YHHPlanarityEnergyCreator ) )
{}

/// clone
EnergyMethodOP
YHHPlanarityEnergy::clone() const
{
	return new YHHPlanarityEnergy;
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

///
void
YHHPlanarityEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const
{
	using numeric::constants::d::degrees_to_radians;
	using numeric::constants::d::pi;

	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ){
			return;
	}
	if ( defines_score_for_rsd(rsd) ) {
		emap[ yhh_planarity ] += 0.5 * (std::cos( pi - 2*rsd.chi(3)*degrees_to_radians)+1);
	}
}


bool
YHHPlanarityEnergy::defines_dof_derivatives( pose::Pose const & ) const
{
	return true;
}

Real
YHHPlanarityEnergy::eval_residue_dof_derivative(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,// min_data,
	id::DOF_ID const &,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const &,// pose,
	ScoreFunction const &,// sfxn,
	EnergyMap const & weights
) const
{
	using numeric::constants::d::degrees_to_radians;
	using numeric::constants::d::pi;
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ){
			return 0.0;
	}

	if ( ! tor_id.valid() ) return 0.0;
	if ( defines_score_for_rsd(rsd) && tor_id.type() == id::CHI && tor_id.torsion() == 3 ) {
		return
			weights[ yhh_planarity ] * std::sin(pi -2*rsd.chi(3)*degrees_to_radians);
	} else {
		return 0.0;
	}
}

/// @brief P_AA_pp_Energy is context independent; indicates that no
/// context graphs are required
void
YHHPlanarityEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}

core::Size
YHHPlanarityEnergy::version() const
{
	return 1; // Initial versioning
}

bool
YHHPlanarityEnergy::defines_score_for_rsd( conformation::Residue const & rsd ) const
{
	return rsd.aa() == chemical::aa_tyr && rsd.type().nchi() == 3 && rsd.type().atom_is_hydrogen( rsd.type().chi_atoms(3)[4] );
}

} // methods
} // scoring
} // core

