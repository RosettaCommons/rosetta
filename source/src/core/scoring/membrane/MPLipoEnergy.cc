// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/membrane/MPLipoEnergy.cc
///
///	@brief		Membrane Lipophibicity Term
///	@details	Whole Structure Energy - Evaluate structure based on derived
///				lipophobicities from input in lips file.
///				Last Modified: 7/6/14
///
///	@author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPLipoEnergy_cc
#define INCLUDED_core_scoring_membrane_MPLipoEnergy_cc

// Unit headers
#include <core/scoring/membrane/MPLipoEnergy.hh> 
#include <core/scoring/membrane/MPLipoEnergyCreator.hh> 

// Project Headers
#include <core/scoring/membrane/MembraneData.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/LipidAccInfo.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "core.scoring.membrane.MPLipoEnergy" );

namespace core {
namespace scoring {
namespace membrane {

using namespace core::scoring;
using namespace core::scoring::methods;

// Creator Methods //////////////////////////////////////////////////

/// @brief Return a fresh instance of the MPLipo energy term
methods::EnergyMethodOP
MPLipoEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
	) const {
	return new MPLipoEnergy;
}

/// @brief Return MPLipo Score type associated with method
ScoreTypes
MPLipoEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( MPLipo );
	return sts;
}

// Constructors /////////////////////////////////////////////////

/// @brief Constructor
MPLipoEnergy::MPLipoEnergy() :
	parent( new MPLipoEnergyCreator ),
	mpdata_( ScoringManager::get_instance()->get_MembraneData() )
{}


/// @brief Clone Energy Method
EnergyMethodOP
MPLipoEnergy::clone() const
{
	return new MPLipoEnergy();
}

// Scoring Methods /////////////////////////////////////////////////

/// @brief Setup Method for Scoring
void
MPLipoEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	// compute interpolated number of neighbors at various distance cutoffs
	pose.update_residue_neighbors();
	mpdata_.compute_centroid_environment( pose );
	
}
	
/// @brief Compute whole structure lipophilicity energy from predicted
/// values in membrane info
void
MPLipoEnergy::finalize_total_energy(
									pose::Pose & pose,
									ScoreFunction const &,
									EnergyMap & emap
									) const {
	
	// Check Structure is a membrane protein
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Error: Cannot use mpframework energy term using a non membrane pose!" );
	}
	
	// Check that the user has included lipid accessibility data
	if ( !include_lips_ ) return; // TODO implement this check
	
	// Obtain a cenlist from the psoe
	CenListInfo const & cenlist = mpdata_.get_cenlist_from_pose( pose );
	
	// Determine total number of tmhs from membrane pose
	TR << "WARNING: MPLipoEnergy: Using the number of spans instead of the number of inserted helices!!!" << std::endl;
	Size num_tmh = pose.conformation().membrane_info()->spanning_topology()->total_spans();
	
	// Initialize MP Lips Score
	Real mplipo_energy( 0 );
	
	// Initialize Component Scores
	Real cen10Buried( 0 );
	Real cen10Exposed( 0 );
	Real cen10Buried_norm( 0 );
	Real cen10Exposed_norm( 0 );

	// Loop through the pose and compute energies
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
				
		// Make Exceptions for specific types of residues
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;
		
		// Otherwise, compute the score for the given residue
		Real B( pose.conformation().membrane_info()->lipid_acc_data()->lipid_exposure()[ i ] );
		Real E( pose.conformation().membrane_info()->lipid_acc_data()->lipid_burial()[ i ] );
		
		if ( B != 0 ) {
			cen10Buried += B * cenlist.fcen10(i);
			cen10Buried_norm += 1;
		}
		if ( E != 0 ) {
			cen10Exposed += E * cenlist.fcen10(i);
			cen10Exposed_norm += 1;
		}
	}
	
	// Compute Lips update
	Real B_mean( 0 );
	Real E_mean( 0 );
	if( cen10Exposed_norm !=0 ) E_mean = cen10Exposed/cen10Exposed_norm;
	if( cen10Buried_norm != 0 ) B_mean = cen10Buried/cen10Buried_norm;
	mplipo_energy += (E_mean-B_mean)*num_tmh;
	
	// Update score in the energy map and finalize whole structure energy
	emap[ MPLipo ] = mplipo_energy;
	mpdata_.finalize( pose );
}

core::Size
MPLipoEnergy::version() const
{
	return 1; // Initial versioning
}
	
} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPLipoEnergy_cc
