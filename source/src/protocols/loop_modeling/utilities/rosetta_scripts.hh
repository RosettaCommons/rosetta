// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_utilities_rosetta_scripts_HH
#define INCLUDED_protocols_loop_modeling_utilities_rosetta_scripts_HH

// Protocol headers
#include <protocols/loop_modeling/LoopMover.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>

// RosettaScripts headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

/// @brief Instantiate a LoopMover from its rosetta-scripts tag.
/// @details If the given tag is not a LoopMover, a nice error message is 
/// thrown.
LoopMoverOP loop_mover_from_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose);

/// @brief Parse the "scorefxn" rosetta-scripts tag for the given mover.
/// @details The given mover must implement a set_score_function() method.
template <class LoopMoverSubclass>
void set_scorefxn_from_tag(
		LoopMoverSubclass & mover,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data) {

	using protocols::rosetta_scripts::parse_score_function;

	// Only call set_score_function() if the tag is actually specified, otherwise 
	// score functions specified in parent movers will be erroneously ignored.

	if (tag->hasOption("scorefxn")) {
		mover.set_score_function(parse_score_function(tag, data, ""));
	}
}

/// @brief Parse the "task_operations" rosetta-scripts tag for the given mover.
/// @details The given mover must implement a set_score_function() method.
template <class LoopMoverSubclass>
void set_task_factory_from_tag(
		LoopMoverSubclass & mover,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data) {

	using protocols::rosetta_scripts::parse_task_operations;

	// Only call set_task_factory() if the tag is actually specified, otherwise 
	// task factories specified in parent movers will be erroneously ignored.

	if (tag->hasOption("task_operations")) {
		mover.set_task_factory(parse_task_operations(tag, data));
	}
}

}
}
}

#endif
