// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobDistributor.cc
/// @brief  JobDistributor class definition
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_JobDistributor_HH
#define INCLUDED_protocols_jd3_JobDistributor_HH

// Unit headers
#include <protocols/jd3/JobDistributor.fwd.hh>

// Package headers
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

namespace protocols {
namespace jd3 {

class JobDistributor : public utility::pointer::ReferenceCount {
public:
	typedef utility::vector1< LarvalJobOP > LarvalJobVector;

public:

	JobDistributor();
	virtual ~JobDistributor();

	/// @brief The main method for executing a protocol.  Derived classes have the option of
	/// overriding the behavior of the base class, or changing the behavior of the base class
	/// version of this function by changing the other virtual methods it calls.
	virtual
	void
	go( JobQueenOP queen );

protected:

	/// @brief Set the job queen -- this must be called by the derived classes if they override
	/// the go() method.
	void set_job_queen( JobQueenOP job_queen );

	/// @brief Store the list of jobs to be run in the current round
	void store_jobs_for_current_round( LarvalJobs const & jobs );

	/// @brief Access to the JobQueen object for derived JobDistributors
	JobQueen & job_queen();

	virtual
	LarvalJobs
	determine_jobs_for_next_round();

	virtual
	bool
	more_jobs_in_current_round();

	virtual
	LarvalJobOP
	select_next_job();

	virtual
	void
	purge_similar_jobs_which_have_bad_inputs( LarvalJobOP job );

	/// @brief Invoke this non-virtual method to iterate across the jobs_for_current_round_ list
	/// stored in the base class; useful if the derived class looks at or asks the
	/// base class to look at the jobs_for_current_round_ list.
	void
	mark_similar_jobs_which_have_bad_inputs_in_job_list( LarvalJobOP job );

	virtual
	void
	process_exception_from_job( LarvalJobOP job, utility::excn::EXCN_Base const & exception );

	virtual
	void
	process_job_result( LarvalJobOP job, JobResultOP result );

	virtual
	void
	note_round_completed();

	virtual
	bool
	another_round_remains();

private:

	JobQueenOP job_queen_;

	LarvalJobVector jobs_for_current_round_;
	core::Size njobs_for_round_;
	core::Size next_job_index_;

};

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd3_JobDistributor_HH
