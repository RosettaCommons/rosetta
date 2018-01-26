// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/PoseResource.hh
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_core_import_pose_PoseResource_HH
#define INCLUDED_core_import_pose_PoseResource_HH

//unit headers
#include <core/import_pose/PoseResource.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

// Basic headers
#include <basic/resource_manager/ResourceLoader.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

//C++ headers
#include <istream>

namespace core {
namespace import_pose {

class PoseResource : public utility::pointer::ReferenceCount
{
public:
	PoseResource( core::pose::PoseCOP pose );
	virtual ~PoseResource();

	void pose( core::pose::PoseCOP pose_to_hold );

	/// @brief The primary way of getting a Pose from a PoseResource:
	/// to request a deep copy of the Pose that it is holding.
	/// This ensures that the Pose is safely copied to prevent thread-unsafe
	/// modification of its internal data.
	///
	/// @details Pose's operator= is not threadsafe in that Poses maintain
	/// pointers to each other to silently communicate between themselves
	/// when one of them is modified, and these connections are
	/// established during an operator= call. The alternative is to call
	/// Pose::detached_copy( Pose const & source ) which this function does.
	/// No thread-unsafe connections are established between Poses copied
	/// with the detached_copy call.
	///
	/// Note that every call to this function creates a new Pose, and
	/// so there might end up with many copies of the Pose in memory. Use
	/// the protocols::moves::PoseFromPoseResourceMover to load a
	/// single Pose from the ResourceManager into the DataMap during
	/// its initialization. Then subsequent movers can read that Pose
	/// from the DataMap during their initialization.
	core::pose::PoseOP pose_deep_copy() const;

private:
	core::pose::PoseCOP pose_;

};

} // namespace import_pose
} // namespace core

#endif //INCLUDED_core_import_pose_PoseResource_HH

