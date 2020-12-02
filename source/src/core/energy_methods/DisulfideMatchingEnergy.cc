// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/DisulfideMatchingEnergy.cc
/// @brief  Centroid Disulfide Energy class implementation
/// @author rvernon@u.washington.edu
/// @date   02/09/10


// Unit headers
#include <core/energy_methods/DisulfideMatchingEnergy.hh>
#include <core/energy_methods/DisulfideMatchingEnergyCreator.hh>

// Package headers
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>
#include <core/scoring/disulfides/DisulfideMatchingEnergyContainer.hh>

// Project headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/Methods.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {


/// @details This must return a fresh instance of the DisulfideMatchingEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
DisulfideMatchingEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< DisulfideMatchingEnergy >( core::scoring::ScoringManager::get_instance()->get_DisulfideMatchingPotential() );
}

core::scoring::ScoreTypes
DisulfideMatchingEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( dslfc_rot );
	sts.push_back( dslfc_trans );
	sts.push_back( dslfc_RT );
	return sts;
}


static basic::Tracer TR( "core.energy_methods.DisulfideMatchingEnergy" );


DisulfideMatchingEnergy::DisulfideMatchingEnergy(
	core::scoring::disulfides::DisulfideMatchingPotential const & potential
) :
	parent( utility::pointer::make_shared< DisulfideMatchingEnergyCreator >() ),
	potential_( potential )
{}

DisulfideMatchingEnergy::~DisulfideMatchingEnergy() = default;

// EnergyMethod Methods:

core::scoring::methods::EnergyMethodOP DisulfideMatchingEnergy::clone() const
{
	return utility::pointer::make_shared< DisulfideMatchingEnergy >( potential_ );
}


void DisulfideMatchingEnergy::setup_for_scoring(
	pose::Pose & pose,
	core::scoring::ScoreFunction const & ) const
{
	using namespace core::scoring::methods;

	if ( pose.energies().long_range_container( disulfide_matching_energy ) == nullptr ) {
		core::scoring::disulfides::DisulfideMatchingEnergyContainerOP dec(
			utility::pointer::make_shared< core::scoring::disulfides::DisulfideMatchingEnergyContainer >( pose )
		);
		pose.energies().set_long_range_container( disulfide_matching_energy, dec );
	} else {
		core::scoring::disulfides::DisulfideMatchingEnergyContainerOP dec(
			utility::pointer::static_pointer_cast< core::scoring::disulfides::DisulfideMatchingEnergyContainer > (
			pose.energies().nonconst_long_range_container( disulfide_matching_energy )
			)
		);
		dec->update( pose );
	}
}

void DisulfideMatchingEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & ) const
{}


// TwoBodyEnergy Methods:

void DisulfideMatchingEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd1.has_variant_type( core::chemical::REPLONLY ) || rsd2.has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}

	Energy match_rot;
	Energy match_trans;
	Energy match_RT;

	//Require cysteines
	//Require Centroid (NOT NO MORE!)
	if ( ( ! rsd1.type().is_sidechain_thiol() && ! rsd1.type().is_disulfide_bonded() ) || (! rsd2.type().is_sidechain_thiol() && ! rsd2.type().is_disulfide_bonded() ) ) return;

	core::scoring::disulfides::DisulfideMatchingEnergyContainerCOP dec(
		utility::pointer::static_pointer_cast< core::scoring::disulfides::DisulfideMatchingEnergyContainer const > (
		pose.energies().long_range_container( core::scoring::methods::disulfide_matching_energy )
		)
	);
	//Require they're bonded
	if ( ! dec->residue_forms_disulfide( rsd1.seqpos() ) ||
			dec->other_neighbor_id( rsd1.seqpos() ) != (Size) rsd2.seqpos() ) {
		return;
	}

	potential_.score_disulfide(
		rsd1, rsd2,
		match_rot,
		match_trans,
		match_RT
	);

	emap[ core::scoring::dslfc_rot ]   += match_rot;
	emap[ core::scoring::dslfc_trans ] += match_trans;
	emap[ core::scoring::dslfc_RT ]    += match_RT;
}


bool DisulfideMatchingEnergy::defines_intrares_energy( core::scoring::EnergyMap const & ) const
{
	return false;
}


void DisulfideMatchingEnergy::eval_intrares_energy(
	conformation::Residue const &,
	pose::Pose const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap &
) const
{}

// LongRangeTwoBodyEnergy methods
core::scoring::methods::LongRangeEnergyType
DisulfideMatchingEnergy::long_range_type() const
{
	return core::scoring::methods::disulfide_matching_energy;
}


bool DisulfideMatchingEnergy::defines_residue_pair_energy(
	pose::Pose const & pose,
	Size rsd1,
	Size rsd2
) const
{
	using namespace core::scoring::methods;
	if ( ! pose.energies().long_range_container( disulfide_matching_energy ) ) return false;

	core::scoring::disulfides::DisulfideMatchingEnergyContainerCOP dec(
		utility::pointer::static_pointer_cast< core::scoring::disulfides::DisulfideMatchingEnergyContainer const > (
		pose.energies().long_range_container( disulfide_matching_energy )
		)
	);
	return dec->disulfide_bonded( rsd1, rsd2 );
}
core::Size
DisulfideMatchingEnergy::version() const
{
	return 1; // Initial versioning
}


} // namespace energy_methods
} // namespace core

