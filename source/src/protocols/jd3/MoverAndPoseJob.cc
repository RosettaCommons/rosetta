// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/MoverAndPoseJob.cc
/// @brief  The class method definitions for MoverAndPoseJob and PoseJobResult
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/MoverAndPoseJob.hh>

// Project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

namespace protocols {
namespace jd3 {

MoverAndPoseJob::MoverAndPoseJob() {}
MoverAndPoseJob::~MoverAndPoseJob() {}

JobResultOP MoverAndPoseJob::run() {
	mover_->apply( *pose_ );
	PoseJobResultOP result = create_job_result();
	result->pose( pose_ );
	finalize_results( result );
	return result;
}

void
MoverAndPoseJob::mover( moves::MoverOP setting ) {
	mover_ = setting;
}

void
MoverAndPoseJob::pose( core::pose::PoseOP setting ) {
	pose_ = setting;
}

moves::MoverOP MoverAndPoseJob::mover() {
	return mover_;
}

core::pose::PoseOP MoverAndPoseJob::pose() {
	return pose_;
}

PoseJobResultOP MoverAndPoseJob::create_job_result() { return PoseJobResultOP( new PoseJobResult ); }
void MoverAndPoseJob::finalize_results( PoseJobResultOP ) {}


PoseJobResult::PoseJobResult() {}
PoseJobResult::~PoseJobResult() {}

JobStatus PoseJobResult::status() const { return jd3_job_status_success; }

core::pose::PoseOP PoseJobResult::pose() { return pose_; }
void PoseJobResult::pose( core::pose::PoseOP setting ) { pose_ = setting; }

} // namespace jd3
} // namespace protocols

