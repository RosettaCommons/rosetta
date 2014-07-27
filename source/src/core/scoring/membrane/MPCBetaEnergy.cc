// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/membrane/MPCbetaEnergy.cc
///
///	@brief		Membrane Environemnt CBeta Energy
///	@details	One Body Term - Score packing density in the membrane. Scores centroids for within
///				6A and 12A radius. Derived from Membrane base potential and uses mpframework data
///				Last Modified: 4/2/14
///
///	@author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPCbetaEnergy_cc
#define INCLUDED_core_scoring_membrane_MPCbetaEnergy_cc

// Unit Headers
#include <core/scoring/membrane/MPCBetaEnergy.hh> 
#include <core/scoring/membrane/MPCBetaEnergyCreator.hh> 

// Project Headers
#include <core/scoring/membrane/MembraneData.hh> 
#include <core/scoring/ScoringManager.hh> 

// Package Headers
#include <core/pose/Pose.hh> 
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/conformation/membrane/MembraneInfo.hh> 
#include <core/conformation/membrane/SpanningTopology.hh>

#include <core/kinematics/Jump.hh>
#include <core/scoring/EnergyMap.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR("core.scoring.membrane.MPCbetaEnergy");

using namespace core::scoring;
using namespace core::scoring::methods;

namespace core {
namespace scoring {
namespace membrane {
	
/// @brief Return a Fresh Instance of the MPCbeta Energy Class
methods::EnergyMethodOP
MPCbetaEnergyCreator::create_energy_method(
										   methods::EnergyMethodOptions const &
										   ) const {
	return new MPCbetaEnergy;
}

/// @brief Return Applicable Score Type for Method (MPCbeta)
ScoreTypes
MPCbetaEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( MPCbeta );
	return sts;
}


/// @brief Default Constructor
MPCbetaEnergy::MPCbetaEnergy() :
	parent( new MPCbetaEnergyCreator ),
	mpdata_( ScoringManager::get_instance()->get_MembraneData() )
{}

/// @brief Clone Method
EnergyMethodOP
MPCbetaEnergy::clone() const {
	return new MPCbetaEnergy;
}

/// @brief Setup Centroid and Mmebrane Potential for Scoring
void
MPCbetaEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	mpdata_.compute_centroid_environment( pose );
}

/// @brief Evlauate CBeta Residue Energy
void
MPCbetaEnergy::residue_energy(
			   conformation::Residue const & rsd,
			   pose::Pose const & pose,
			   EnergyMap & emap
				) const {
	
	// Check Structure is a membrane protein
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message("Error: Cannot use mpframework energy term using a non membrane pose!");
	}
	
	// Initialize Scores
	Real membrane_cb_score( 0.0 );
	
	// Skip Virtual Residues
	if ( rsd.aa() == core::chemical::aa_vrt ) {
		emap[ MPCbeta ] = membrane_cb_score;
		return;
	}
			
	// Initialize Info for Scoring
	CenListInfo const & cenlist = mpdata_.get_cenlist_from_pose( pose );
	core::Size const seqpos = rsd.seqpos();
	
	// Determine total number of tmhs from membrane pose
	Size num_tmh = pose.conformation().membrane_info()->spanning_topology()->total_spans();
	
	// Compute CBeta score at the position
	membrane_cb_score = compute_mpcbeta_score( cenlist, seqpos, num_tmh );
	membrane_cb_score = 2.667 * membrane_cb_score * 0.3;
	emap[ MPCbeta ] += membrane_cb_score;
	
}

/// @brief Finalize Total Energy
void
MPCbetaEnergy::finalize_total_energy(
					  pose::Pose & pose,
					  ScoreFunction const &,
					  EnergyMap &
					  ) const
{
	mpdata_.finalize( pose );
}

/// @brief Specify Version for Energy Method
core::Size
MPCbetaEnergy::version() const { return 2; }

/// @brief Evaluate CBeta Energy Term = Actual Evaluate
core::Real
MPCbetaEnergy::compute_mpcbeta_score(
									 CenListInfo const & cenlist,
									 core::Size const seqpos,
									 core::Size const num_tmh
									 ) const
{
	// Initialize Cbeta Score
	core::Real mp_cbeta_score( 0.0 );
	core::Real mp_cbeta_score6( 0.0 );
	core::Real mp_cbeta_score12( 0.0 );

	// Compute fcen6 and fcen12 from database
	Real const fcen6( cenlist.fcen6(seqpos) );
	Real const fcen12( cenlist.fcen12(seqpos) );
	
	// interp1 rounds down to nearest (non-negative) integer.
	// note cen6 is always at least 1.0
	int const interp1 = static_cast< int >( fcen6 );
	int const interp3 = static_cast< int >( fcen12 );
	
	// lower bound
	Real const interp2 = fcen6-interp1;
	Real const interp4 = fcen12-interp3;
	
	TR.Debug << "WARNING: MPCBetaEnergy does not currently depend on the number of inserted helices!!!" << std::endl;
	
	// Evaluate score for 6A radius from num tmhs
	if ( num_tmh <= 2 ) {
		
		mp_cbeta_score6 =
		( 1.0-interp2 ) * mpdata_.mem_cbeta_2TM_den6()( interp1 )+
		interp2         * mpdata_.mem_cbeta_2TM_den6()( interp1+1 );
		
	} else if ( num_tmh <= 4 ) {
		
		mp_cbeta_score6 =
		(1.0-interp2) * mpdata_.mem_cbeta_4TM_den6()( interp1 )+
		interp2       * mpdata_.mem_cbeta_4TM_den6()( interp1+1 );
		
	} else {
		
		mp_cbeta_score6 =
		(1.0-interp2) * mpdata_.mem_cbeta_den6()( interp1 )+
		interp2       * mpdata_.mem_cbeta_den6()( interp1+1 );
	}
	
	// Evaluate Score for 12A radius
	mp_cbeta_score12 =
	(1.0-interp4) * mpdata_.mem_cbeta_den12()( interp3 )+
	interp4       * mpdata_.mem_cbeta_den12()( interp3+1 );
	
	mp_cbeta_score = ( mp_cbeta_score6 + mp_cbeta_score12 );
	return mp_cbeta_score;
}

	
} // membrane
} // scoring
} // core


#endif // INCLUDED_core_scoring_membrane_MPCbetaEnergy_cc
