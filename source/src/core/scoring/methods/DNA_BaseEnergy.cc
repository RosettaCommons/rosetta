// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/DNA_BaseEnergy.cc
/// @brief
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/DNA_BaseEnergy.hh>
#include <core/scoring/methods/DNA_BaseEnergyCreator.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/dna/DNA_BasePotential.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/dna/setup.hh> // set_base_partner

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


// ObjexxFCL headers

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the DNA_BaseEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
DNA_BaseEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new DNA_BaseEnergy );
}

ScoreTypes
DNA_BaseEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( dna_bp );
	sts.push_back( dna_bs );
	return sts;
}


using namespace dna; //////////////////// NOTE NOTE NOTE

/// the atom through which the knowledge based potential applies a force

std::string const dna_deriv_atom( " C5 " );

DNA_BaseEnergy::DNA_BaseEnergy() :
	parent( methods::EnergyMethodCreatorOP( new DNA_BaseEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_DNA_BasePotential() )
{}


/// clone
EnergyMethodOP
DNA_BaseEnergy::clone() const
{
	return EnergyMethodOP( new DNA_BaseEnergy() );
}

/// are these really necessary??????????? move to scheme that doesnt depend on nbr calcn
///
void
DNA_BaseEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	set_base_partner( pose );
	pose.update_residue_neighbors();
}


void
DNA_BaseEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	set_base_partner( pose );
	pose.update_residue_neighbors();
}

void
DNA_BaseEnergy::setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const
{
	set_base_partner( pose );
	pose.update_residue_neighbors();
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// ordered!!!! requires pos1<pos2
inline
bool
count_pair_bs(
	Size const pos1,
	Size const pos2,
	BasePartner const & partner
)
{
	return ( pos2 == pos1 + 1 && partner[pos1] && partner[pos2] && partner[pos2] == partner[pos1]-1 &&
		partner[pos1] != pos2 );
}


/// same as dna::retrieve_base_partner_from_pose
inline
BasePartner const &
retrieve_base_partner_from_pose_inline( pose::Pose const & pose )
{
	//using core::pose::datacache::CacheableDataType::BASE_PARTNER;
	debug_assert( pose.data().has( core::pose::datacache::CacheableDataType::BASE_PARTNER ) );
	debug_assert( dynamic_cast< BasePartner const *>( &( pose.data().get( core::pose::datacache::CacheableDataType::BASE_PARTNER ))));
	return ( static_cast< BasePartner const &>(    pose.data().get( core::pose::datacache::CacheableDataType::BASE_PARTNER )));
}


void
DNA_BaseEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	if ( !rsd1.is_DNA() || !rsd2.is_DNA() ) return;

	Real bp_score( 0.0 ), bs_score( 0.0 );

	// retrieve DNA basepair info from pose
	BasePartner const & base_partner( retrieve_base_partner_from_pose_inline( pose ) );

	// base step score:
	Size const pos1( rsd1.seqpos() ), pos2( rsd2.seqpos() );
	if ( count_pair_bs( pos1, pos2, base_partner ) ) {
		bs_score += potential_.base_step_score( rsd1, rsd2 );
	}
	if ( count_pair_bs( pos2, pos1, base_partner ) ) {
		bs_score += potential_.base_step_score( rsd2, rsd1 );
	}


	// base pair score
	if ( pos2 == base_partner[ pos1 ] ) {
		if ( pos1 < pos2 ) {
			bp_score += potential_.base_pair_score( rsd1, rsd2 );
		} else {
			bp_score += potential_.base_pair_score( rsd2, rsd1 );
		}
	}

	emap[ dna_bs ] += bs_score;
	emap[ dna_bp ] += bp_score;
	//std::cout << "DNA_BaseEnergy " << rsd1.seqpos() << " " << rsd2.seqpos() << " " << bs_score << " " << bp_score << std::endl;
}


void
DNA_BaseEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	Size const  pos1( atom_id.rsd() );
	Size const atom1( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( pos1 ) );
	if ( ( !rsd1.is_DNA() ) || ( rsd1.atom_name( atom1 ) != dna_deriv_atom ) ) return;

	// retrieve DNA basepair info from pose
	BasePartner const & base_partner( retrieve_base_partner_from_pose_inline( pose ) );

	///////////////////
	// base step derivs
	if ( weights[ dna_bs ] != 0.0 ) {
		// to next residue:
		if ( pos1 < pose.total_residue() ) {
			Size const pos2( pos1+1 );
			conformation::Residue const & rsd2( pose.residue( pos2 ) );
			if ( rsd2.is_DNA() && count_pair_bs( pos1, pos2, base_partner ) ) {
				potential_.eval_base_step_derivative( rsd1, rsd2, F1, F2,        weights[ dna_bs ] );
			}
		}

		// to previous residue:
		if ( pos1 > 1 ) {
			Size const pos2( pos1-1 );
			conformation::Residue const & rsd2( pose.residue( pos2 ) );
			if ( rsd2.is_DNA() && count_pair_bs( pos2, pos1, base_partner ) ) {
				potential_.eval_base_step_derivative( rsd2, rsd1, F1, F2, -1.0 * weights[ dna_bs ] );
			}
		}
	}

	///////////////////
	// base pair derivs
	if ( weights[ dna_bp ] != 0.0 && base_partner[ pos1 ] ) {
		Size const pos2( base_partner[ pos1 ] );
		if ( pos1 < pos2 ) {
			potential_.eval_base_pair_derivative( rsd1, pose.residue( pos2 ), F1, F2,        weights[ dna_bp ] );
		} else {
			potential_.eval_base_pair_derivative( pose.residue( pos2 ), rsd1, F1, F2, -1.0 * weights[ dna_bp ] );
		}
	}
}


/// @brief DNA_BaseEnergy distance cutoff
Distance
DNA_BaseEnergy::atomic_interaction_cutoff() const
{
	return 5.5; // -- temporary hack to allow us to use the standard neighbor array
}

/// @brief DNA_BaseEnergy
void
DNA_BaseEnergy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}
core::Size
DNA_BaseEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core
