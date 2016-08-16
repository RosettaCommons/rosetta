// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author James Thompson

// libRosetta headers

#ifndef INCLUDED_core_import_pose_pose_stream_ExtendedPoseInputStream_HH
#define INCLUDED_core_import_pose_pose_stream_ExtendedPoseInputStream_HH

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/import_pose/pose_stream/PoseInputStream.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace core {
namespace import_pose {
namespace pose_stream {

class ExtendedPoseInputStream : public PoseInputStream {

public:
	ExtendedPoseInputStream( std::string sequence, Size ntimes )
	: seq_( sequence ), ntimes_( ntimes ), current_n_( 1 ) {}
	~ExtendedPoseInputStream() {}

	virtual bool has_another_pose();

	virtual void reset();

	virtual void fill_pose(
		core::pose::Pose & pose,
		core::chemical::ResidueTypeSet const & residue_set
	);
	virtual void fill_pose( core::pose::Pose& );

private:
	std::string seq_;
	Size ntimes_, current_n_;
}; // ExtendedPoseInputStream

} // pose_stream
} // import_pose
} // core

#endif
