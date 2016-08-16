// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/copy_dofs/CopyDofStepWiseSampler.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/copy_dofs/CopyDofStepWiseSampler.hh>
#include <protocols/simple_moves/CopyDofMover.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.copy_dofs.CopyDofStepWiseSampler" );

namespace protocols {
namespace stepwise {
namespace sampler {
namespace copy_dofs {

//Constructor
CopyDofStepWiseSampler::CopyDofStepWiseSampler(utility::vector1< core::pose::PoseOP > const & pose_list,
	std::map< Size, Size > const & res_map,
	core::pose::Pose const & starting_pose )
{
	for ( Size n = 1; n <= pose_list.size(); n++ ) {
		simple_moves::CopyDofMover copy_dof_mover( *pose_list[n], res_map );
		core::pose::PoseOP pose_copy = starting_pose.clone();
		copy_dof_mover.apply( *pose_copy );
		pose_list_.push_back( pose_copy );
		copy_dof_movers_.push_back( simple_moves::CopyDofMoverOP( new simple_moves::CopyDofMover( *pose_copy, res_map ) ) );
	}
}

//Constructor
CopyDofStepWiseSampler::CopyDofStepWiseSampler(utility::vector1< core::pose::PoseOP > const & pose_list,
	std::map< Size, Size > const & res_map ):
	pose_list_( pose_list )
{
	for ( Size n = 1; n <= pose_list.size(); n++ ) copy_dof_movers_.push_back( simple_moves::CopyDofMoverOP( new simple_moves::CopyDofMover( *pose_list[n], res_map ) ) );
}

///////////////////////////////////////////////////
void
CopyDofStepWiseSampler::apply( core::pose::Pose & pose, core::Size const i ){
	copy_dof_movers_[ i ]->apply( pose );
}

} //copy_dofs
} //sampler
} //stepwise
} //protocols
