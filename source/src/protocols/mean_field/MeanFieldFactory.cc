// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/InteractionGraphFactory.cc
/// @brief  Interation graph factory class definition
/// @author Aliza Rubenstein (aliza.rubenstein@gmail.com)

// Unit headers
#include <protocols/mean_field/MeanFieldFactory.hh>

// Package headers
#include <protocols/mean_field/MeanField.hh>
#include <protocols/mean_field/DesignMeanField.hh>
#include <protocols/mean_field/FlexBBMeanField.hh>
#include <protocols/mean_field/FlexBBDesignMeanField.hh>


#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>


#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace mean_field {

/// @details chooses type of MeanField based on whether designable task and how many poses are present
MeanFieldOP
MeanFieldFactory::create_mean_field(
	core::Size option,
	core::pose::PoseOPs & poses,
	utility::vector1 < core::pack::task::PackerTaskOP > tasks,
	core::scoring::ScoreFunctionOP scfxn,
	core::Real lambda_mem,
	core::Real conv_diff,
	core::Real temp,
	core::Real threshold
)
{

	if ( ! tasks[ 1 ]->design_any() ) {
		if ( poses.size() == 1 ) {
			MeanFieldOP mf( new MeanField ( option, poses, tasks, scfxn, lambda_mem, conv_diff, temp, threshold ) );
			return mf;
		} else {
			FlexBBMeanFieldOP fbmf( new FlexBBMeanField ( option, poses, tasks, scfxn, lambda_mem, conv_diff, temp, threshold ) );
			return fbmf;
		}
	} else {
		if ( poses.size() == 1 ) {
			DesignMeanFieldOP dmf( new DesignMeanField ( option, poses, tasks, scfxn, lambda_mem, conv_diff, temp, threshold ) );
			return dmf;
		} else {
			FlexBBDesignMeanFieldOP fbdmf( new FlexBBDesignMeanField ( option, poses, tasks, scfxn, lambda_mem, conv_diff, temp, threshold ) );
			return fbdmf;
		}
	}

}

}
}

