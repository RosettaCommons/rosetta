// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/job_distributors/MultiThreadedJobDistributor.cc
/// @brief  MultiThreadedJobDistributor class definition
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_job_distributors_MultiThreadedJobDistributor_HH
#define INCLUDED_protocols_jd3_job_distributors_MultiThreadedJobDistributor_HH


#ifdef MULTI_THREADED
#ifdef CXX11

// Unit headers
#include <protocols/jd3/job_distributors/MultiThreadedJobDistributor.fwd.hh>

// Package headers
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/JobQueen.fwd.hh>
#include <protocols/jd3/Job.fwd.hh>
#include <protocols/jd3/JobResult.fwd.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.fwd.hh>

// C++ headers
#include <atomic>

namespace protocols {
namespace jd3 {
namespace job_distributors {

class MultiThreadedJobDistributor : public JobDistributor {
public:
	typedef utility::vector1< LarvalJobOP > LarvalJobVector;

public:

	MultiThreadedJobDistributor();
	virtual ~MultiThreadedJobDistributor();

	/// @brief The main method for executing a protocol.  Derived classes have the option of
	/// overriding the behavior of the base class, or changing the behavior of the base class
	/// version of this function by changing the other virtual methods it calls.
	virtual
	void
	go( JobQueenOP queen );

protected:

	/// @brief Poll the running jobs, and if there are any that have completed,
	/// process their output, and then return true.
	bool
	capacity_for_another_thread();

	bool
	unfinished_jobs_remain();

	void
	process_finished_job( core::Size ii );

	void
	start_job( LarvalJobOP, JobOP );

private:

	core::Size nthreads_;
	core::Size count_since_last_job_finished_;
	utility::vector1< JobRunnerOP > jobs_running_in_each_thread_;

};

class JobRunner {
public:
	JobRunner( LarvalJobOP, JobOP );
	void run();
	bool exited_w_exception() const;
	std::string const & exception_message() const;
	LarvalJobOP larval_job();
	JobOP mature_job();
	JobResultOP job_result();
	bool complete();
private:
	LarvalJobOP larval_job_;
	JobOP       mature_job_;
	JobResultOP job_result_;

	std::string exception_message_;
	bool        exited_w_exception_;

	std::atomic< bool > complete_;
};

} // job_distributors
} // jd3
} // protocols

#endif // MULTI_THREADED

#endif // INCLUDED_protocols_jd3_job_distributors_MultiThreadedJobDistributor_HH
