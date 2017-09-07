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

// Unit headers
#include <protocols/jd3/job_distributors/MultiThreadedJobDistributor.fwd.hh>

// Package headers
#include <protocols/jd3/CompletedJobOutput.hh>
#include <protocols/jd3/Job.fwd.hh>
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/JobDigraph.fwd.hh>
#include <protocols/jd3/JobQueen.fwd.hh>
#include <protocols/jd3/JobResult.fwd.hh>
#include <protocols/jd3/JobSummary.fwd.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>
#include <protocols/jd3/job_distributors/JobExtractor.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.fwd.hh>

// C++ headers
#include <atomic>
#include <map>

// CTPL headers
#include <CTPL/ctpl_stl.fwd.h>

namespace protocols {
namespace jd3 {
namespace job_distributors {

class MultiThreadedJobDistributor : public JobDistributor {
public:
	typedef utility::vector1< LarvalJobOP > LarvalJobVector;

	typedef std::list< core::Size > SizeList;
	typedef std::map< JobResultID, std::pair< LarvalJobOP, JobResultOP > > JobResultMap;

	typedef utility::pointer::shared_ptr< ctpl::thread_pool > ThreadPoolOP;

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

	bool
	prepare_next_job( LarvalJobOP larval_job, core::Size attempt_count );

	bool
	unfinished_jobs_remain();

	void
	process_completed_job( JobRunnerOP runner );

	void potentially_output_some_job_results();
	void potentially_discard_some_job_results();

	void
	start_job( LarvalJobOP, JobOP );

private:

	JobQueenOP job_queen_;
	JobDigraphOP job_dag_;
	JobExtractorOP job_extractor_;

	// The big old map that stores all the JobResults that are generated
	// over the course of execution
	JobResultMap job_results_;

	core::Size nthreads_;
	std::list< JobRunnerOP > jobs_running_;

	Size default_retry_limit_;

	ThreadPoolOP thread_pool_;

};

class JobRunner {
public:
	JobRunner( LarvalJobOP, JobOP, core::Size attempt_count, core::Size retry_limit );

	/// @brief The main function for this unit of work.
	void run( int thread_index );

	bool exited_w_exception() const;
	std::string const & exception_message() const;
	LarvalJobOP larval_job() const;
	JobOP mature_job() const;
	core::Size attempt_count() const;

	CompletedJobOutput const & job_output() const;
	bool complete() const;
	int running_thread() const;
private:
	core::Size attempt_count_;
	core::Size retry_limit_;
	int  running_thread_;

	LarvalJobOP larval_job_;
	JobOP       mature_job_;

	CompletedJobOutput job_result_;

	std::string exception_message_;
	bool        exited_w_exception_;

	std::atomic< bool > complete_;
};

} // job_distributors
} // jd3
} // protocols

#endif // MULTI_THREADED

#endif // INCLUDED_protocols_jd3_job_distributors_MultiThreadedJobDistributor_HH
