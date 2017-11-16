// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/legacy/screener/SimplePoseSelection.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/legacy/screener/SimplePoseSelection.hh>
#include <protocols/stepwise/modeler/align/StepWiseLegacyClusterer.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.legacy.screener.SimplePoseSelection" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace screener {

//Constructor
SimplePoseSelection::SimplePoseSelection( pose::Pose const & pose,
	utility::vector1< Size > const & moving_res_list,
	modeler::options::StepWiseModelerOptionsCOP options,
	bool const full_optimize ):
	pose_( pose ),
	moving_res_list_( moving_res_list ),
	options_( options ),
	full_optimize_( full_optimize )
{}

//Destructor
SimplePoseSelection::~SimplePoseSelection()
{}

bool
SimplePoseSelection::check_screen() {
	pose_list_.push_back( pose_.clone() );
	return true;
}

void
SimplePoseSelection::finalize() {
	if ( pose_list_.size() == 0 ) return;
	modeler::align::StepWiseLegacyClusterer stepwise_clusterer( pose_list_, moving_res_list_,
		options_, full_optimize_ /* force align*/ );
	stepwise_clusterer.cluster();
	pose_list_ = stepwise_clusterer.get_pose_list();
}

} //screener
} //legacy
} //stepwise
} //protocols
