// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/methods/DNA_ReferenceEnergy.cc
/// @brief  dna scoring
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/DNA_ReferenceEnergy.hh>
#include <core/scoring/methods/DNA_ReferenceEnergyCreator.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/conformation/Residue.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh> // HACK

// Utility headers


// C++


namespace core {
namespace scoring {
namespace methods {

using namespace dna; //////////////////// NOTE NOTE NOTE
static basic::Tracer TR("core.scoring.methods.DNA_ReferenceEnergy");



methods::EnergyMethodOP
DNA_ReferenceEnergyCreator::create_energy_method(
	EnergyMethodOptions const & options
) const {
	return EnergyMethodOP( new DNA_ReferenceEnergy( options ) );
}

ScoreTypes
DNA_ReferenceEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( dna_ref );
	return sts;
}

DNA_ReferenceEnergy::DNA_ReferenceEnergy( EnergyMethodOptions const & options ):
	parent( EnergyMethodCreatorOP( new DNA_ReferenceEnergyCreator ) )
{
	if ( !options.has_method_weights( dna_ref ) ) utility_exit_with_message( "dna_ref requires method weights!" );
	// order is aa,ac,ag,at,ca,cc,cg,ct,...
	utility::vector1< Real > wts( options.method_weights( dna_ref ) );
	base_step_reference_energies_.resize(4);
	runtime_assert( wts.size() == 16 );
	for ( Size i=1; i<= 4; ++i ) {
		base_step_reference_energies_[i].resize(4);
		for ( Size j=1; j<= 4; ++j ) {
			base_step_reference_energies_[i][j] = wts.front();
			wts.erase( wts.begin() );
		}
	}
	runtime_assert( wts.empty() );
}

/// copy c-tor
DNA_ReferenceEnergy::DNA_ReferenceEnergy( DNA_ReferenceEnergy const & src ):
	parent( EnergyMethodCreatorOP( new DNA_ReferenceEnergyCreator ) ),
	base_step_reference_energies_( src.base_step_reference_energies_ )
{
}


/// clone
EnergyMethodOP
DNA_ReferenceEnergy::clone() const
{
	return EnergyMethodOP( new DNA_ReferenceEnergy( *this ) );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


/// same as dna::retrieve_base_partner_from_pose
inline
BasePartner const &
retrieve_base_partner_from_pose_inline( pose::Pose const & pose )
{
	assert( pose.data().has( core::pose::datacache::CacheableDataType::BASE_PARTNER ) );
	assert( dynamic_cast< BasePartner const *>( &( pose.data().get( core::pose::datacache::CacheableDataType::BASE_PARTNER ))));
	return ( static_cast< BasePartner const &>(    pose.data().get( core::pose::datacache::CacheableDataType::BASE_PARTNER )));
}


///
/// 04/29/11 === change to allow rsd1 and rsd2 in either sequence order
///
void
DNA_ReferenceEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	if ( !rsd1.is_DNA() || !rsd2.is_DNA() ) return;

	Size const pos1( rsd1.seqpos() ), pos2( rsd2.seqpos() );
	if ( pos2 == pos1 + 1 && !rsd1.is_upper_terminus() &&
			count_pair_bs( pos1, pos2, retrieve_base_partner_from_pose_inline( pose ) ) ) {
		emap[ dna_ref ] += base_step_energy( rsd1.aa(), rsd2.aa() );
	} else if ( pos1 == pos2+1 && !rsd2.is_upper_terminus() &&
			count_pair_bs( pos2, pos1, retrieve_base_partner_from_pose_inline( pose ) ) ) {
		emap[ dna_ref ] += base_step_energy( rsd2.aa(), rsd1.aa() );
	}
}

/// @brief DNA_ReferenceEnergy distance cutoff
Distance
DNA_ReferenceEnergy::atomic_interaction_cutoff() const
{
	return 5.5; // -- temporary hack to allow us to use the standard neighbor array
}

/// @brief DNA_ReferenceEnergy
void
DNA_ReferenceEnergy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}



} // methods
} // scoring
} // core
