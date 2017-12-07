// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/import_pose/pose_stream/MetaPoseInputStream.hh
/// @brief
/// @author James Thompson

// libRosetta headers

#ifndef INCLUDED_core_import_pose_pose_stream_MetaPoseInputStream_HH
#define INCLUDED_core_import_pose_pose_stream_MetaPoseInputStream_HH

#include <core/types.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/import_pose/pose_stream/PoseInputStream.hh>

#include <utility/vector1.fwd.hh>


namespace core {
namespace import_pose {
namespace pose_stream {

class MetaPoseInputStream : public PoseInputStream {

public:
	MetaPoseInputStream() : current_index_( 1 ) {}

	~MetaPoseInputStream() {}

	void add_pose_input_stream( PoseInputStreamOP input );

	utility::vector1< PoseInputStreamOP > get_input_streams();

	virtual bool has_another_pose();

	virtual void reset();

	virtual void fill_pose(
		core::pose::Pose & pose,
		core::chemical::ResidueTypeSet const & residue_set
	);
	virtual void fill_pose( core::pose::Pose& );

private:
	core::Size current_index_;
	utility::vector1< PoseInputStreamOP > input_streams_;
	utility::vector1< PoseInputStreamOP >::iterator current_input_stream_;
}; // MetaPoseInputStream

} // pose_stream
} // import_pose
} // core

#endif
