// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/legacy/screener/SimplePoseSelection.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/legacy/screener/SimplePoseSelection.hh>
#include <protocols/stepwise/sampling/align/StepWiseLegacyClusterer.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.legacy.screener.SimplePoseSelection" );

namespace protocols {
namespace stepwise {
namespace legacy {
namespace screener {

	//Constructor
	SimplePoseSelection::SimplePoseSelection( pose::Pose const & pose,
																						utility::vector1< Size > const & moving_res_list,
																						sampling::modeler_options::StepWiseModelerOptionsCOP options,
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
		sampling::align::StepWiseLegacyClusterer stepwise_clusterer( pose_list_,	moving_res_list_,
																														 options_, full_optimize_ /* force align*/ );
		stepwise_clusterer.cluster();
		pose_list_ = stepwise_clusterer.get_pose_list();
	}

} //screener
} //legacy
} //stepwise
} //protocols
