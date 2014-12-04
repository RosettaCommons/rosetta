// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopMoverGroup_HH
#define INCLUDED_protocols_loop_modeling_LoopMoverGroup_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/utilities/LoopMoverGroup.fwd.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocols headers
#include <protocols/filters/Filter.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

/// @brief Manage a collection of loop-sampling algorithms.
///
/// @details This class is meant to make a group of LoopMover objects behave
/// like a single LoopMover.  Methods like add_mover() and add_filter() are
/// also provided to add movers to the group.

class LoopMoverGroup : public LoopMover {

public:

	/// @brief Default constructor.
	LoopMoverGroup();

	/// @brief Return the name of this mover.
	string get_name() const { return "LoopMoverGroup"; }

	/// @brief Add the names of all the algorithms invoked by this loop mover to 
	/// the given list.  Indentation is used to represent hierarchy.
	void get_children_names(
			utility::vector1<std::string> & names, std::string indent="") const;

	/// @brief Add a mover to this group.
	LoopMoverOP add_mover(LoopMoverOP task);

	/// @brief Add a filter to this group.
	LoopMoverOP add_filter(protocols::filters::FilterOP filter);

	/// @brief Add a nested group to this group.
	LoopMoverGroupOP add_mover_group();

	/// @brief Remove all movers and filters from this group.
	void clear();

	/// @brief Indicate that the current set of movers in this group is meant as
	/// some sort of default, and should be cleared when a new mover is added.
	void mark_as_default();

	/// @brief Return true if no movers are contained in this group.
	bool empty() const;

protected:
	/// @brief Apply all the movers added to this group.
	/// @details This method will immediately return false if any of the child
	/// movers reports failure.
	bool do_apply(Pose & pose);

private:
	bool is_default_;

};

}
}
}

#endif
