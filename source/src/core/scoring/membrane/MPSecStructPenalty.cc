// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/membrane/MPSecStructPenalty.cc
///
///	@brief		Penalty for Non Helix Secondary Structures in the Membrane
///				based on PsiPred prediction
///	@details	One Body Term - evaluate secondary structure residue by
///				residue and sum the penalty.
///				Last Modified: 4/2/14
///
///	@author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPSecStructPenalty_cc
#define INCLUDED_core_scoring_membrane_MPSecStructPenalty_cc

// Unit Headers
#include <core/scoring/membrane/MPSecStructPenalty.hh>
#include <core/scoring/membrane/MPSecStructPenaltyCreator.hh>

// Project Headers
#include <core/scoring/membrane/MembraneData.hh>
#include <core/scoring/ScoringManager.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/scoring/EnergyMap.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>

using namespace core::scoring;
using namespace core::scoring::methods;

namespace core {
namespace scoring {
namespace membrane {
	
/// Creator Methods ////////////////////

/// @brief Return a Fresh Instance of the Energy Method
methods::EnergyMethodOP
MPSecStructPenaltyCreator::create_energy_method(
										 methods::EnergyMethodOptions const &
										 ) const {
	return new MPSecStructPenalty;
}

/// @brief Return Appropriate Score Type Name
ScoreTypes
MPSecStructPenaltyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( MPSecStuc );
	return sts;
}

/// @brief Default Constructor
MPSecStructPenalty::MPSecStructPenalty() :
	parent( new MPSecStructPenaltyCreator ),
	mpdata_( ScoringManager::get_instance()->get_MembraneData() )
{}

/// @brief Clone an Energy Method
EnergyMethodOP
MPSecStructPenalty::clone() const
{
	return new MPSecStructPenalty;
}

/// Scoring Methods //////////////////////

/// @brief Setup Method for Scoring
void
MPSecStructPenalty::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	mpdata_.compute_centroid_environment( pose );
}

/// @brief Setup Energy Method for Derivatives
void
MPSecStructPenalty::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sf) const
{
	setup_for_scoring( pose, sf );
}

/// @brief Score Residue Energy Method (Calls internal evaluate env method)
void
MPSecStructPenalty::residue_energy(
							conformation::Residue const & rsd,
							pose::Pose const & pose,
							EnergyMap & emap
							) const
{
	
	// Check Structure is a membrane protein
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message("Error: Cannot use mpframework energy term using a non membrane pose!");
	}
	
	// wouldn't it be cool if this were layer based? do we have statistis on this?
	// it could be a new term 
	

	// Set the new score in the energies map
	emap[ MPSecStruc ] += env_score;
	
} // residue_energy

/// @brief Finalize Total Energy in the Poose
void
MPSecStructPenalty::finalize_total_energy(
								   pose::Pose & pose,
								   ScoreFunction const &,
								   EnergyMap &
								   ) const
{
	mpdata_.finalize( pose );
}

core::Size
MPSecStructPenalty::version() const
{
	return 2;
}
		
} // membrane
} // scoring
} // core


#endif // INCLUDED_core_scoring_membrane_MPSecStructPenalty_cc
