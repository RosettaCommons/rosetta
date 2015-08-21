// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/scoring_util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_scoring_util_HH
#define INCLUDED_protocols_stepwise_modeler_scoring_util_HH

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.fwd.hh>

namespace protocols {
namespace stepwise {
namespace modeler {

core::scoring::ScoreFunctionOP
get_minimize_scorefxn( core::pose::Pose const & pose,
	core::scoring::ScoreFunctionCOP scorefxn,
	options::StepWiseModelerOptionsCOP options );


core::scoring::ScoreFunctionCOP
initialize_sample_scorefxn( core::scoring::ScoreFunctionCOP scorefxn,
	core::pose::Pose const & pose,
	options::StepWiseModelerOptionsCOP options );

core::scoring::ScoreFunctionCOP
initialize_pack_scorefxn( core::scoring::ScoreFunctionCOP sample_scorefxn, core::pose::Pose const & pose );

core::scoring::ScoreFunctionCOP
initialize_o2prime_pack_scorefxn( core::scoring::ScoreFunctionCOP const & scorefxn );

} //modeler
} //stepwise
} //protocols

#endif
