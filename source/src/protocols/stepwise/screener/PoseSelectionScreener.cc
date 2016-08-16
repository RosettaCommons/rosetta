// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/PoseSelectionScreener.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/PoseSelectionScreener.hh>
#include <protocols/stepwise/modeler/align/StepWiseClusterer.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.screener.PoseSelectionScreener" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {

//Constructor
PoseSelectionScreener::PoseSelectionScreener( pose::Pose & pose,
	scoring::ScoreFunctionCOP scorefxn,
	modeler::align::StepWiseClustererOP stepwise_clusterer ):
	pose_( pose ),
	scorefxn_( scorefxn ),
	stepwise_clusterer_( stepwise_clusterer )
{
}

//Destructor
PoseSelectionScreener::~PoseSelectionScreener()
{}

bool
PoseSelectionScreener::check_screen(){
	( *scorefxn_ )( pose_ );
	stepwise_clusterer_->apply( pose_ ); // super-simple. sticks into clusterer.
	return true;
}

} //screener
} //stepwise
} //protocols
