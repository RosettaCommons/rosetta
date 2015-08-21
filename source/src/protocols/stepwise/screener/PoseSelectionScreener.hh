// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/PoseSelectionScreener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_PoseSelectionScreener_HH
#define INCLUDED_protocols_stepwise_screener_PoseSelectionScreener_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/modeler/align/StepWiseClusterer.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace screener {

class PoseSelectionScreener: public StepWiseScreener {

public:

	//constructor
	PoseSelectionScreener( core::pose::Pose & pose,
		core::scoring::ScoreFunctionCOP scorefxn,
		modeler::align::StepWiseClustererOP stepwise_clusterer );

	//destructor
	~PoseSelectionScreener();

public:

	std::string
	name() const { return "PoseSelectionScreener"; }

	StepWiseScreenerType
	type() const { return POSE_SELECTION; }

	bool
	check_screen();

private:

	core::pose::Pose & pose_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	modeler::align::StepWiseClustererOP stepwise_clusterer_;

};

} //screener
} //stepwise
} //protocols

#endif
