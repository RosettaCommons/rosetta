// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/MoverAndPoseJob.hh
/// @brief  The definition of class protocols::jd3::MoverAndPoseJob, which applies a
/// Mover to a Pose in its run() method.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_MoverAndPoseJob_HH
#define INCLUDED_protocols_jd3_MoverAndPoseJob_HH

// Unit headers
#include <protocols/jd3/MoverAndPoseJob.fwd.hh>

// Package headers
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/JobResult.hh>

// Project headers
#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace jd3 {

class MoverAndPoseJob : public Job
{
public:

	MoverAndPoseJob();
	~MoverAndPoseJob() override;

	
	JobResultOP run() override;

	void mover( moves::MoverOP setting );
	void pose( core::pose::PoseOP setting );

	moves::MoverOP mover();
	core::pose::PoseOP pose();

protected:
	/// @brief Factory method so derived classes can return their own result classes,
	/// which themselves should derive from PoseJobResult.
	virtual PoseJobResultOP create_job_result();

	/// @brief Method that allows derived classes to tuck data into the result object
	/// as they see fit.
	virtual void finalize_results( PoseJobResultOP result );

private:
	moves::MoverOP mover_;
	core::pose::PoseOP pose_;

};

class PoseJobResult : public JobResult
{
public:
	PoseJobResult();
	~PoseJobResult() override;

	JobStatus status() const override;

	core::pose::PoseOP pose();
	void pose( core::pose::PoseOP setting );
private:

	core::pose::PoseOP pose_;

};

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_Job_HH
