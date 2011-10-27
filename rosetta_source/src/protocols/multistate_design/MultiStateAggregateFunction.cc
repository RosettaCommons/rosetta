// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MultiStateAggregateFunction.tmpl.hh
/// @brief
/// @author Colin A. Smith

// Unit headers
#include <protocols/multistate_design/MultiStateAggregateFunction.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh> // REQUIRED FOR WINDOWS

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
// Auto-header: duplicate removed #include <protocols/multistate_design/MultiStateAggregateFunction.hh>
#include <protocols/multistate_design/MultiStateFitnessFunction.hh>
#include <protocols/multistate_design/SingleState.fwd.hh>
#include <protocols/multistate_design/SingleState.hh>
#include <protocols/multistate_design/SingleStateFitnessFunction.fwd.hh>
#include <protocols/multistate_design/SingleStateFitnessFunction.hh>
#include <protocols/viewer/viewers.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.fwd.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <iostream>


namespace protocols {
namespace multistate_design {

core::Real
MultiStateAggregateFunction::evaluate(
	utility::vector1<core::Real> const & single_state_fitnesses,
	MultiStateFitnessFunction & fitness_function
) const
{
	utility::vector1<SingleStateCOP> single_states(fitness_function.const_states());
	runtime_assert(single_state_fitnesses.size() == single_states.size());

	core::Real aggregate_fitness = 0;

	for (core::Size i = 1; i <= single_state_fitnesses.size(); ++i) {
		aggregate_fitness += (single_states[i]->is_positive_state() ? 1.0 : -1.0) * single_state_fitnesses[i];
	}

	return aggregate_fitness;
}

} // namespace multistate_design
} // namespace protocols


