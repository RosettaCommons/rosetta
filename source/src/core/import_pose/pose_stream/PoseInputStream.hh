// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PoseInputStream.hh
/// @brief
/// @author James Thompson


#ifndef INCLUDED_core_import_pose_pose_stream_PoseInputStream_HH
#define INCLUDED_core_import_pose_pose_stream_PoseInputStream_HH

#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace import_pose {
namespace pose_stream {

class PoseInputStream : public utility::pointer::ReferenceCount {
public:
	PoseInputStream() {}
	virtual ~PoseInputStream() {}

	virtual bool has_another_pose() = 0;

	virtual void reset() = 0;

	virtual void fill_pose(
		core::pose::Pose &,
		core::chemical::ResidueTypeSet const &,
		bool const metapatches = true
	) = 0;

	virtual void fill_pose( core::pose::Pose&, bool const metapatches = true ) = 0;

	virtual utility::vector1< core::pose::PoseOP > get_all_poses(
		core::chemical::ResidueTypeSet const & residue_set
	);

	/// @brief Perform common operations on a Pose dependent on command-line
	/// options before finishing fill_pose. This includes optimizing hydrogens,
	/// adding constraints, and finding disulfides. This should be called in the
	/// fill_pose method of derived classes.
	void preprocess_pose( core::pose::Pose & pose );

}; // class PoseInputStream

} // pose_stream
} // import_pose
} // core

#endif
