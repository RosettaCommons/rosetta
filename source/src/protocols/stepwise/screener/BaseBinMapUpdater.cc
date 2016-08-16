// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/BaseBinMapUpdater.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/BaseBinMapUpdater.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueList.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.screener.BaseBinMapUpdater" );

using namespace protocols::stepwise::sampler;
using namespace protocols::stepwise::sampler::rigid_body;

namespace protocols {
namespace stepwise {
namespace screener {

//Constructor
BaseBinMapUpdater::BaseBinMapUpdater( BaseBinMap & base_bin_map ):
	base_bin_map_( base_bin_map )
{}

//Destructor
BaseBinMapUpdater::~BaseBinMapUpdater()
{}

// diagnostics.
void
BaseBinMapUpdater::get_update( sampler::StepWiseSamplerBaseOP sampler ){

	if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_LIST ) {
		RigidBodyStepWiseSamplerWithResidueList & rigid_body_rotamer_with_copy_dofs = *( static_cast< RigidBodyStepWiseSamplerWithResidueList * >( sampler.get() ) );
		update_base_bin_map( rigid_body_rotamer_with_copy_dofs.get_rigid_body_values() );
		return;
	}
	if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ) {
		RigidBodyStepWiseSamplerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyStepWiseSamplerWithResidueAlternatives * >( sampler.get() ) );
		update_base_bin_map( rigid_body_rotamer_with_residue_alternatives.get_rigid_body_values() );
		return;
	}

	runtime_assert( sampler->type() == RIGID_BODY );
	RigidBodyStepWiseSampler & rigid_body_rotamer = *( static_cast< RigidBodyStepWiseSampler * >( sampler.get() ) );
	update_base_bin_map( rigid_body_rotamer.get_values() );
}


/////////////////////////////////////////////////////////////////////////////////////
// diagnostics
void
BaseBinMapUpdater::update_base_bin_map( BaseBin const & base_bin ){
	std::map< BaseBin, int, compare_base_bin > ::const_iterator it = base_bin_map_.find( base_bin );
	if ( it == base_bin_map_.end() ) base_bin_map_[base_bin] = 0;
	base_bin_map_[base_bin] ++;
}

/////////////////////////////////////////////////////////////////////////////////////
// diagnostics
void
BaseBinMapUpdater::update_base_bin_map( utility::vector1< Real > const & rigid_body_values ){
	BaseBin base_bin;
	base_bin.centroid_x  = (int)(rigid_body_values[6]);
	base_bin.centroid_y  = (int)(rigid_body_values[5]);
	base_bin.centroid_z  = (int)(rigid_body_values[4]);
	base_bin.euler_alpha = (int)(rigid_body_values[3]);
	base_bin.euler_z     = (int)(rigid_body_values[2]);
	base_bin.euler_gamma = (int)(rigid_body_values[1]);
	update_base_bin_map( base_bin );
}


} //screener
} //stepwise
} //protocols
