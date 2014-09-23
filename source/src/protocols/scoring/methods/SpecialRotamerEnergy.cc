// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/scoring/methods/SpecialRotamerEnergy.hh
/// @brief  Adds a bonus to any rotamer that is flagged
/// @author sthyme, sthyme@gmail.com, Feb 2010

// Unit headers
#include <protocols/scoring/methods/SpecialRotamerEnergy.hh>
#include <protocols/scoring/methods/SpecialRotamerEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/EnergyMap.hh>

// Project headers
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>



namespace protocols {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the SpecialRotamerEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
SpecialRotamerEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return core::scoring::methods::EnergyMethodOP( new SpecialRotamerEnergy );
}

core::scoring::ScoreTypes
SpecialRotamerEnergyCreator::score_types_for_method() const {
	core::scoring::ScoreTypes sts;
	sts.push_back( core::scoring::special_rot );
	return sts;
}


/// ctor
SpecialRotamerEnergy::SpecialRotamerEnergy() :
	parent( core::scoring::methods::EnergyMethodCreatorOP( new SpecialRotamerEnergyCreator ) )
{}

/// clone
core::scoring::methods::EnergyMethodOP
SpecialRotamerEnergy::clone() const
{
	return core::scoring::methods::EnergyMethodOP( new SpecialRotamerEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

///
void
SpecialRotamerEnergy::residue_energy(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & ,
	core::scoring::EnergyMap & emap
) const
{
	if ( rsd.has_variant_type( core::chemical::SPECIAL_ROT ) ) {
		emap[ core::scoring::special_rot ] = 1.0;
	}
}


///
core::Real
SpecialRotamerEnergy::eval_dof_derivative(
	core::id::DOF_ID const &,// dof_id,
	core::id::TorsionID const &, //  tor_id
	core::pose::Pose const &, // pose
	core::scoring::ScoreFunction const &, //sfxn
	core::scoring::EnergyMap const & // weights
) const
{
	return 0.0;
}

/// @brief SpecialRotamerEnergy is context independent; indicates that no
/// context graphs are required
void
SpecialRotamerEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}
core::Size
SpecialRotamerEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // protocols
