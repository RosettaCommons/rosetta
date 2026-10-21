// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rosetta_scripts/BatchedMultiplePoseMover.hh
/// @brief Base class for a mover that processes all poses produced by the previous mover at once.
///   Before apply(), ParsedProtocol drains the previous mover's get_additional_output()
///   and passes the pose handed to apply() followed by every drained pose to set_batched_poses().
///   Subclasses read them with batched_poses() inside apply().
/// @author Federico Olivieri (federico.alberto.olivieri@ki.si)

#ifndef INCLUDED_protocols_rosetta_scripts_BatchedMultiplePoseMover_hh
#define INCLUDED_protocols_rosetta_scripts_BatchedMultiplePoseMover_hh

#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace rosetta_scripts {

class BatchedMultiplePoseMover;
typedef utility::pointer::shared_ptr< BatchedMultiplePoseMover > BatchedMultiplePoseMoverOP;
typedef utility::pointer::shared_ptr< BatchedMultiplePoseMover const > BatchedMultiplePoseMoverCOP;

class BatchedMultiplePoseMover : public protocols::moves::Mover {

public:
	~BatchedMultiplePoseMover() override = default;

	/// @brief Used by RosettaScripts to set the poses to process: the input pose first, then the previous mover's additional outputs
	void set_batched_poses( utility::vector1< core::pose::PoseOP > const & poses ) { batched_poses_ = poses; }

	utility::vector1< core::pose::PoseOP > const & batched_poses() const { return batched_poses_; }

private:
	utility::vector1< core::pose::PoseOP > batched_poses_;

};

} // rosetta_scripts
} // protocols

#endif
