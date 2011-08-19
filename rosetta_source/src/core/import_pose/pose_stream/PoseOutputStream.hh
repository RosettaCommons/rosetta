// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/import_pose/pose_stream/PoseOutputStream.hh
/// @brief
/// @author James Thompson


#ifndef core_io_pose_stream_PoseOutputStream_HH
#define core_io_pose_stream_PoseOutputStream_HH

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>

namespace core {
namespace import_pose {
namespace pose_stream {

class PoseOutputStream : public utility::pointer::ReferenceCount {

public:

	PoseOutputStream()  {}
	~PoseOutputStream() {}

	virtual void write_pose( core::pose::Pose & pose ) = 0;

}; // class PoseOutputStream

//typedef utility::pointer::owning_ptr< PoseOutputStream > PoseOutputStreamOP;

} // pose_stream
} // import_pose
} // core

#endif
