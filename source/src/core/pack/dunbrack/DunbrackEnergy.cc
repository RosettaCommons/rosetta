// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/dunbrack/DunbrackEnergy.cc
/// @brief  Dunbrack energy method implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/dunbrack/DunbrackEnergy.hh>
#include <core/pack/dunbrack/DunbrackEnergyCreator.hh>

// Package Headers
#include <core/id/PartialAtomID.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>

#include <core/scoring/ScoreType.hh>

// Project headers
#include <core/chemical/rotamers/NCAARotamerLibrarySpecification.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <core/id/TorsionID.hh>

// Utility headers
#include <numeric/conversions.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace dunbrack {

using namespace scoring;
using namespace scoring::methods;

/// @details This must return a fresh instance of the DunbrackEnergy class,
/// never an instance already in use
scoring::methods::EnergyMethodOP
DunbrackEnergyCreator::create_energy_method(
	scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< DunbrackEnergy >();
}

scoring::ScoreTypes
DunbrackEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_dun );
	sts.push_back( fa_dun_dev );
	sts.push_back( fa_dun_rot );
	sts.push_back( fa_dun_semi );
	return sts;
}

/// ctor
DunbrackEnergy::DunbrackEnergy() :
	parent( utility::pointer::make_shared< DunbrackEnergyCreator >() )
{}

DunbrackEnergy::~DunbrackEnergy() = default;

/// clone
EnergyMethodOP
DunbrackEnergy::clone() const
{
	return utility::pointer::make_shared< DunbrackEnergy >();
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

/// @details Allocate the scratch space object on the stack to
/// alieviate thread-safety concerns.  Scratch does not use new.
void
DunbrackEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	scoring::EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) )  return;

	if ( rsd.is_virtual_residue() ) return;

	//Returns the equivalent L-amino acid library if a D-amino acid is provided
	pack::rotamers::SingleResidueRotamerLibraryCOP rotlib = rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rsd.type() );

	if ( ! rotlib || rsd.has_variant_type( core::chemical::SC_BRANCH_POINT ) ) return;

	rotamers::TorsionEnergy tenergy;
	emap[ fa_dun ] += rotlib->rotamer_energy( rsd, pose, tenergy );
	emap[ fa_dun_rot ] += tenergy.rot;
	emap[ fa_dun_semi ] += tenergy.semi;
	emap[ fa_dun_dev  ] += tenergy.dev;
}

bool DunbrackEnergy::defines_dof_derivatives( pose::Pose const & ) const { return true; }

utility::vector1< id::PartialAtomID >
DunbrackEnergy::atoms_with_dof_derivatives(
	conformation::Residue const & res,
	pose::Pose const & pose
) const
{
	utility::vector1< id::PartialAtomID > retlist;
	pack::rotamers::SingleResidueRotamerLibraryCOP rotlib =
		rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( res.type() );
	if ( ! rotlib || res.has_variant_type( core::chemical::SC_BRANCH_POINT ) )  { return retlist; }

	std::set< id::PartialAtomID > atoms = rotlib->atoms_w_dof_derivatives(res, pose);

	retlist.resize(atoms.size());
	std::copy(atoms.begin(), atoms.end(), retlist.begin());
	return retlist;
}


Real
DunbrackEnergy::eval_residue_dof_derivative(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & ,//min_data,
	id::DOF_ID const & ,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	scoring::ScoreFunction const & ,//sfxn,
	scoring::EnergyMap const & weights
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) return 0.0;

	if ( ! tor_id.valid() ) return 0.0;

	pack::rotamers::SingleResidueRotamerLibraryCOP rotlib =
		rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rsd.type() );

	if ( ! rotlib || rsd.has_variant_type( core::chemical::SC_BRANCH_POINT ) )  { return 0.0; }

	rotamers::TorsionEnergy tderiv;
	rotlib->rotamer_energy_deriv( rsd, pose, tor_id, tderiv );

	return numeric::conversions::degrees( weights[ fa_dun ] * tderiv.tot + weights[ fa_dun_dev ] * tderiv.dev + weights[ fa_dun_rot ] * tderiv.rot + weights[ fa_dun_semi ] * tderiv.semi);
}

/// @brief DunbrackEnergy is context independent; indicates that no context graphs are required
void
DunbrackEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
) const
{}

core::Size
DunbrackEnergy::version() const
{
	return 1; // Initial versioning
}


} // dunbrack
} // pack
} // core

