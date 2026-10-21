// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rosetta_scripts/AdditionalPosesMover.hh
/// @brief This mover is an abstract base class for movers which
///   accept multiple poses from a previous mover.
///   (This is for integration with other processing machinery.)
///
///   It's a mover which *processes* additional Poses,
///   versus one which can create addtional Poses (which all movers could do if they want)
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_moves_AdditionalPosesMover_hh
#define INCLUDED_protocols_moves_AdditionalPosesMover_hh

// Package headers
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <string>

namespace protocols {
namespace moves {

class AdditionalPosesMover;
typedef utility::pointer::shared_ptr< AdditionalPosesMover > AdditionalPosesMoverOP;
typedef utility::pointer::shared_ptr< AdditionalPosesMover const > AdditionalPosesMoverCOP;

class AdditionalPosesMover : public protocols::moves::Mover {

public:
	AdditionalPosesMover() = default;
	AdditionalPosesMover(std::string const & name = ""): Mover(name) {}

	// Abstract class -- don't need to (re)implement all of the Mover-necessary functions

	/// @brief Used by RosettaScripts to set the previous mover to pull poses from
	void set_previous_mover( protocols::moves::MoverOP const m );

	/// @brief Test if the previous mover is set
	bool previous_mover_set() const;

	/// @brief The name of the previous mover (for logging purposes)
	/// (Robust to previous mover not being set.)
	std::string previous_mover_name() const;

	/// @brief Get the next pose from the previous mover
	/// May return nullptr if there aren't any to get or if previous_mover isn't set
	/// @details Don't forget the pose passed to apply() as well.
	core::pose::PoseOP get_next_pose() const;

	/// @brief Get the next batch of N poses from the additional output of the previous mover
	/// Pass zero if you want all available poses - Warning! Be careful of previous movers who lazily produce infinite output.
	/// May not return the full N poses if there aren't any to get.
	/// (Will return an empty vector with warning if previous_mover isn't set)
	/// @details Don't forget the pose passed to apply() as well.
	utility::vector1< core::pose::PoseOP > get_next_pose_batch(core::Size N) const;

private:
	protocols::moves::MoverOP previous_mover_ = nullptr;
};

} //moves
} //protocols

#endif
