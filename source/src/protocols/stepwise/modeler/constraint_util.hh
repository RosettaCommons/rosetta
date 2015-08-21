// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/constraint_util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_constraint_util_HH
#define INCLUDED_protocols_stepwise_modeler_constraint_util_HH

#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace modeler {

core::scoring::constraints::ConstraintSetOP
constraint_set_slice( core::scoring::constraints::ConstraintSetOP & cst_set,
	utility::vector1< core::Size > const & slice_res,
	core::pose::Pose const & pose,
	core::pose::Pose const & full_pose );

void
check_scorefxn_has_constraint_terms_if_pose_has_constraints( core::pose::Pose const & pose,
	core::scoring::ScoreFunctionOP & scorefxn );


} //modeler
} //stepwise
} //protocols

#endif
