// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/Scorer.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/Scorer.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <utility>

static basic::Tracer TR( "protocols.stepwise.screener.Scorer" );

namespace protocols {
namespace stepwise {
namespace screener {

//Constructor
Scorer::Scorer():
	pose_( *( new core::pose::Pose ) )
{}

//Constructor
Scorer::Scorer( core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scorefxn ):
	pose_( pose ),
	scorefxn_(std::move( scorefxn ))
{
}

//Destructor
Scorer::~Scorer() = default;

////////////////////////////////////////////////////////
bool
Scorer::check_screen(){
	( *scorefxn_ )( pose_ );
	return true;
}

} //screener
} //stepwise
} //protocols
