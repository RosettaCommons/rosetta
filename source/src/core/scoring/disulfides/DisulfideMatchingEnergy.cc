// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/DisulfideMatchingEnergy.cc
/// @brief  Centroid Disulfide Energy class implementation
/// @author rvernon@u.washington.edu
/// @date   02/09/10


// Unit headers
#include <core/scoring/disulfides/DisulfideMatchingEnergy.hh>
#include <core/scoring/disulfides/DisulfideMatchingEnergyCreator.hh>

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
namespace scoring {
namespace disulfides {


/// @details This must return a fresh instance of the DisulfideMatchingEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
DisulfideMatchingEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new DisulfideMatchingEnergy( ScoringManager::get_instance()->get_DisulfideMatchingPotential() ) );
}

ScoreTypes
DisulfideMatchingEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( dslfc_rot );
	sts.push_back( dslfc_trans );
	sts.push_back( dslfc_RT );
	return sts;
}


static thread_local basic::Tracer TR( "core.scoring.disulfides.DisulfideMatchingEnergy" );



DisulfideMatchingEnergy::DisulfideMatchingEnergy(
	DisulfideMatchingPotential const & potential
) :
	parent( methods::EnergyMethodCreatorOP( methods::EnergyMethodCreatorOP( new DisulfideMatchingEnergyCreator ) ) ),
	potential_( potential )
{}

DisulfideMatchingEnergy::~DisulfideMatchingEnergy() {}

// EnergyMethod Methods:

methods::EnergyMethodOP DisulfideMatchingEnergy::clone() const
{
	return methods::EnergyMethodOP( new DisulfideMatchingEnergy( potential_ ) );
}


void DisulfideMatchingEnergy::setup_for_scoring(
		pose::Pose & pose,
		ScoreFunction const & ) const
{
	using namespace methods;

	if ( pose.energies().long_range_container( disulfide_matching_energy ) == 0 ) {
		DisulfideMatchingEnergyContainerOP dec( new DisulfideMatchingEnergyContainer( pose ) );
		pose.energies().set_long_range_container( disulfide_matching_energy, dec );
	} else {
		DisulfideMatchingEnergyContainerOP dec = DisulfideMatchingEnergyContainerOP (
				utility::pointer::static_pointer_cast< core::scoring::disulfides::DisulfideMatchingEnergyContainer > ( pose.energies().nonconst_long_range_container( disulfide_matching_energy ) ));
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
		ScoreFunction const &,
		EnergyMap & emap
		) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd1.has_variant_type( core::chemical::REPLONLY ) || rsd2.has_variant_type( core::chemical::REPLONLY ) ){
		return;
	}

	Energy match_rot;
	Energy match_trans;
	Energy match_RT;

	//Require cysteines
	if ( rsd1.aa() != chemical::aa_cys || rsd2.aa() != chemical::aa_cys ) return;
	//Require Centroid (NOT NO MORE!)
	//if (rsd1.residue_type_set().name() != chemical::CENTROID ||
	//		rsd2.residue_type_set().name() != chemical::CENTROID )
	//	return;

	DisulfideMatchingEnergyContainerCOP dec = DisulfideMatchingEnergyContainerCOP (
			utility::pointer::static_pointer_cast< core::scoring::disulfides::DisulfideMatchingEnergyContainer const > ( pose.energies().long_range_container( methods::disulfide_matching_energy ) ));
	//Require they're bonded
	if ( ! dec->residue_forms_disulfide( rsd1.seqpos() ) ||
			dec->other_neighbor_id( rsd1.seqpos() ) != (Size) rsd2.seqpos() ){
		return;
	}

	potential_.score_disulfide(
			rsd1, rsd2,
			match_rot,
			match_trans,
			match_RT
			);

	emap[ dslfc_rot ]   += match_rot;
	emap[ dslfc_trans ] += match_trans;
	emap[ dslfc_RT ]    += match_RT;
}



bool DisulfideMatchingEnergy::defines_intrares_energy( EnergyMap const & ) const
{
	return false;
}


void DisulfideMatchingEnergy::eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
		) const
{}

// LongRangeTwoBodyEnergy methods
methods::LongRangeEnergyType
DisulfideMatchingEnergy::long_range_type() const
{
	return methods::disulfide_matching_energy;
}


bool DisulfideMatchingEnergy::defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
		) const
{
	using namespace methods;
	if ( ! pose.energies().long_range_container( disulfide_matching_energy )) return false;

	DisulfideMatchingEnergyContainerCOP dec = DisulfideMatchingEnergyContainerCOP (
			utility::pointer::static_pointer_cast< core::scoring::disulfides::DisulfideMatchingEnergyContainer const > ( pose.energies().long_range_container( disulfide_matching_energy ) ));
	return dec->disulfide_bonded( res1, res2 );
}
core::Size
DisulfideMatchingEnergy::version() const
{
	return 1; // Initial versioning
}


} // namespace disulfides
} // namespace scoring
} // namespace core

