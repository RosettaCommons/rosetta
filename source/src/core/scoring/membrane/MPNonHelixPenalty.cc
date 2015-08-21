// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPNonHelixPenalty.cc
///
/// @brief  Membrane Protein Non helix in Mmebrane Penalty
/// @details Whole structure energy - penalty for helices not in the membrane?
///    and uses mpframework data
///    Last Modified: 3/31/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPNonHelixPenalty_cc
#define INCLUDED_core_scoring_membrane_MPNonHelixPenalty_cc

// Unit Headers
#include <core/scoring/membrane/MPNonHelixPenalty.hh>
#include <core/scoring/membrane/MPNonHelixPenaltyCreator.hh>

// Project Headers
#include <core/scoring/membrane/MembraneData.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/conformation/membrane/SpanningTopology.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/types.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Utility Headers
#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static thread_local basic::Tracer TR( "core.scoring.membrane.MPNonHelix" );

namespace core {
namespace scoring {
namespace membrane {

/// Creator Methods ///////////////////////

/// @brief Return a Fresh Instance of the Energy Method
methods::EnergyMethodOP
MPNonHelixPenaltyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MPNonHelixPenalty );
}

/// @brief Log Score Types
ScoreTypes
MPNonHelixPenaltyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( MPNonHelix );
	return sts;
}


/// @brief Default Constructor
MPNonHelixPenalty::MPNonHelixPenalty() :
	parent( EnergyMethodCreatorOP( new MPNonHelixPenaltyCreator ) ),
	mpdata_( ScoringManager::get_instance()->get_MembraneData() )
{}

/// @brief Clone Method
EnergyMethodOP
MPNonHelixPenalty::clone() const
{
	return EnergyMethodOP( new MPNonHelixPenalty() );
}

/// Scoring Methods /////////////////////////

/// @brief Compute termini penalty per-residue
void
MPNonHelixPenalty::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const {

	using namespace core::scoring::membrane;
	using namespace core::conformation::membrane;
	using namespace core::conformation::symmetry;

	// Check Structure is a membrane protein
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message("Error: Cannot use mpframework energy term using a non membrane pose!");
	}

	// Initialize Penalty Score
	Real non_helix_pen(0);

	// Get Topology from the pose
	SpanningTopologyOP topology = pose.conformation().membrane_info()->spanning_topology();

	// Skip Membrane Residues
	core::Size nres = pose.total_residue()-1;
	if ( rsd.seqpos() > nres ) return;

	// Skip Cases
	if ( rsd.seqpos() == 0 ) return;
	if ( rsd.aa() == core::chemical::aa_vrt ) return;
	// if ( rsd.chain() == topology.size() ) return; // WHAT IS THIS???

	// Grab appropriate resnum and chain num for topology
	core::Size chain = rsd.chain();
	core::Size resnum = rsd.seqpos() - ( pose.conformation().chain_begin(chain) - 1);

	// Compute Info for Penalty
	core::Real z_position = pose.conformation().membrane_info()->residue_z_position( rsd.seqpos() );
	bool tmregion = topology->in_span(resnum);
	char secstruc = pose.conformation().secstruct( rsd.seqpos() );

	// Compute Penalty (using initial z_position, will need to figure out the rsdSeq changes since its only for symm)
	non_helix_pen = compute_nonhelix_penalty( tmregion, secstruc, z_position ) * 10;

	emap[ MPNonHelix ] += non_helix_pen;
}

void
MPNonHelixPenalty::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap &
) const {
	mpdata_.finalize( pose );
}

/// @brief Compute Non Helix in membrane Penalty
core::Real
MPNonHelixPenalty::compute_nonhelix_penalty( bool tmregion, char secstruc, core::Real z_position ) const {

	// Evaluate the penalty
	if ( tmregion && (secstruc != 'H') ) {
		if ( z_position > -12.0 && z_position < 12.0 ) {
			return 1;
		}
	}

	return 0;
}


} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPNonHelixPenalty_cc
