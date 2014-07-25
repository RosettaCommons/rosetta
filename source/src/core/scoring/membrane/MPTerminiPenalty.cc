// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/membrane/MPTerminiPenalty.cc
///
///	@brief		Membrane Protein Termini Penalty
///	@details	Whole structure energy - penalty for residues on the wrong side of the membrane?
///				nd uses mpframework data
///				Last Modified: 3/31/14
///
///	@author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPTerminiPenalty_cc
#define INCLUDED_core_scoring_membrane_MPTerminiPenalty_cc

// Unit Headers
#include <core/scoring/membrane/MPTerminiPenalty.hh> 
#include <core/scoring/membrane/MPTerminiPenaltyCreator.hh> 

// Project Headers
#include <core/scoring/membrane/MembraneData.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>

#include <core/types.hh>

// Utility Headers
#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR("core.scoring.membrane.MPTerminiPenalty");

namespace core {
namespace scoring {
namespace membrane {
	
/// Creator Methods ///////////////////////

/// @brief Return a Fresh Instance of the Energy Method
methods::EnergyMethodOP
MPTerminiPenaltyCreator::create_energy_method(
											  methods::EnergyMethodOptions const &
											  ) const {
	return new MPTerminiPenalty;
}

/// @brief Return Relevant Score Types
ScoreTypes
MPTerminiPenaltyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( MPTermini );
	return sts;
}


/// @brief Default Constructor
MPTerminiPenalty::MPTerminiPenalty() :
	parent( new MPTerminiPenaltyCreator ),
	mpdata_( ScoringManager::get_instance()->get_MembraneData() )
{}

/// @brief Clone Method
EnergyMethodOP
MPTerminiPenalty::clone() const
{
	return new MPTerminiPenalty();
}

/// Scoring Methods /////////////////////////

/// @brief Compute termini penalty per-residue
void
MPTerminiPenalty::residue_energy(
								 conformation::Residue const & rsd,
								 pose::Pose const & pose,
								 EnergyMap & emap
								 ) const {
	
	using namespace core::scoring::membrane;
	
	// Check Structure is a membrane protein
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message("Error: Cannot use mpframework energy term using a non membrane pose!");
	}
	
	// Initialize Termnini Penalty
	Real termini_pen(0);
	
	// Skip Membrane Residues
	core::Size nres = pose.total_residue()-1;
	if ( rsd.seqpos() > nres ) return;

	// Skip Termini Residues
	if ( !rsd.is_terminus() ) return;
	if ( rsd.seqpos() == 0 ) return;
	if ( rsd.aa() == core::chemical::aa_vrt ) return;
	
	// Compute my actual penalty
	core::Real z_position = pose.conformation().residue_z_position( rsd.seqpos() );
	termini_pen = compute_termini_penalty( z_position ) * 50;

	// Add to Score Map
	emap[ MPTermini ] += termini_pen;
}
	
/// @brief Finalize Whole Structure Energy
void
MPTerminiPenalty::finalize_total_energy(
										pose::Pose & pose,
										ScoreFunction const &,
										EnergyMap &
										) const {
	mpdata_.finalize( pose );
}

/// @brief Compute Termini Penalty
/// @brief      Evaluate Termini Penalty
/// @details    Penalty for termini positions contradictory to membrane
///             spanning
core::Real
MPTerminiPenalty::compute_termini_penalty( core::Real z_position ) const {
	
	if ( z_position > -12 && z_position < 12 ) return 1;
	else return 0;
}

	
} // membrane
} // scoring
} // core


#endif // INCLUDED_core_scoring_membrane_MPTerminiPenalty_cc
