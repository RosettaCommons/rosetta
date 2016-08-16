// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/MultiThreadingJob.hh
/// @author James Thompson

#ifndef INCLUDED_protocols_jd2_MultiThreadingJob_hh
#define INCLUDED_protocols_jd2_MultiThreadingJob_hh

//unit headers
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/MultiThreadingJob.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

//C++ headers
#include <string>

#include <protocols/jd2/InnerMultiThreadingJob.fwd.hh>


namespace protocols {
namespace jd2 {

class MultiThreadingJob : public protocols::jd2::Job
{
	MultiThreadingJob(InnerMultiThreadingJobOP inner_job, core::Size nstruct_index);

public:

	InnerMultiThreadingJobOP multi_threading_inner_job();

	virtual ~MultiThreadingJob();

private:
	InnerMultiThreadingJobOP inner_job_;
}; // Job

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_Job_HH
