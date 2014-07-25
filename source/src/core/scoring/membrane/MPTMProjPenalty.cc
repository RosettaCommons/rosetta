// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/membrane/MPTMProjPenalty.cc
///
///	@brief		Membrane Protein TM Proj Penalty
///	@details	Whole structure energy - Penalty for unreasonable tm-helix length compared to predicted
///				helix length (from topology) and uses mpframework data
///				Last Modified: 4/3/14
///
///	@author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPTMProjPenalty_cc
#define INCLUDED_core_scoring_membrane_MPTMProjPenalty_cc

// Unit Headers
#include <core/scoring/membrane/MPTMProjPenalty.hh>
#include <core/scoring/membrane/MPTMProjPenaltyCreator.hh>

// Project Headers
#include <core/scoring/membrane/MembraneData.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/conformation/membrane/MembraneInfo.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/membrane/Span.hh>

#include <core/types.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Utility Headers
#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace scoring {
namespace membrane {
	
/// Creator Methods ///////////////////////

/// @brief Return a Fresh Instance of the Energy Method
methods::EnergyMethodOP
MPTMProjPenaltyCreator::create_energy_method(
											   methods::EnergyMethodOptions const &
											   ) const {
	return new MPTMProjPenalty;
}

/// @brief Return score type MPTMProj
ScoreTypes
MPTMProjPenaltyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( MPTMProj );
	return sts;
}


/// @brief Default Constructor
MPTMProjPenalty::MPTMProjPenalty() :
	parent( new MPTMProjPenaltyCreator ),
	mpdata_( ScoringManager::get_instance()->get_MembraneData() )
{}

/// @brief Clone Method
EnergyMethodOP
MPTMProjPenalty::clone() const
{
	return new MPTMProjPenalty();
}

/// Scoring Methods /////////////////////////

void
MPTMProjPenalty::finalize_total_energy(
										 pose::Pose & pose,
										 ScoreFunction const &,
										 EnergyMap & emap
										 ) const {
	
	using namespace core::conformation::membrane;
	using namespace core::scoring::membrane;
	
	// Initialize TM Projection
	core::Real tm_proj( 0.0 );
	
    // shorthand
    ConformationOP conf = &pose.conformation();
    MembraneInfoOP meminfo = conf->membrane();
    
	// Initialize WHole Pose Data
	core::Vector center = conf->membrane_center();
	core::Vector normal = conf->membrane_normal();

	// Get Topology from the pose
    SpanningTopologyOP topology = meminfo->spanning_topology();
	
	// Read through spanning topology
	for ( Size j = 1; j <= topology->total_spans(); ++j ) {
		
		// Get the center and normal z position
		Real const & start_z_pos = pose.conformation().residue_z_position( topology->span(j)->start() );
		Real const & end_z_pos = pose.conformation().residue_z_position( topology->span(j)->end() );
		
		// Compute seuqnece distance between the two residues
		core::Real dist = topology->span(j)->end()-topology->span(j)->start()+1;
		
		// Compute TM Pcoj penalty for the given helix
		tm_proj += compute_tmproj_penalty( start_z_pos, end_z_pos, dist );
		
	}
	
	// Finalize penaltty
	tm_proj*=50;
	emap[ MPTMProj ] = tm_proj;
	
	// FInalize whole structure energy
	mpdata_.finalize( pose );
}
	
/// @brief Compute Penalty for Length of Helix
core::Real
MPTMProjPenalty::compute_tmproj_penalty(
										core::Real start_z_pos,
										core::Real end_z_pos,
										core::Real dist
										) const {
	
	// Initialize TM projection penalty
	Real tm_proj( 0.0 );
	
	// Calculate helix length
	Real tm_length = std::abs( start_z_pos - end_z_pos );
	
	// Calculate ratio
	Real ratio = tm_length / dist;
	
	// Evaluate penalty
	if(tm_length<15) { tm_proj++; }
	if(ratio<1 || ratio > 1.5) { tm_proj++; }
	
	return tm_proj;
}
	
	
	
} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPTMProjPenalty_cc
