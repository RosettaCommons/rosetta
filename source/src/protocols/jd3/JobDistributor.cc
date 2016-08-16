// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/JobDistributor.cc
/// @brief  JobDistributor class method definitions
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit headers
#include <protocols/jd3/JobDistributor.hh>

// Package headers
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/JobQueen.hh>
#include <protocols/jd3/JobResult.hh>
#include <protocols/jd3/LarvalJob.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.hh>

// C++ headers
#include <string>
#include <ctime>


namespace protocols {
namespace jd3 {

static THREAD_LOCAL basic::Tracer TR( "protocols.jd3.JobDistributor" );

JobDistributor::JobDistributor() :
	njobs_for_round_( 0 ),
	next_job_index_( 1 )
{}

JobDistributor::~JobDistributor() {}

void
JobDistributor::go( JobQueenOP queen ) {
	set_job_queen( queen );

	// something like this should happen -- or be in a function called from this one
	// so that the XSD can be written out without any work being performed.
	// if ( basic::options::option[ basic::options::OptionKeys::jd3::output_xsd ] ) {
	//   TR << queen->job_definition_xsd() << std::endl;
	//   return;
	// }

	do {
		store_jobs_for_current_round( determine_jobs_for_next_round() );
		while ( more_jobs_in_current_round() ) {

			// select the next job to run
			LarvalJobOP larval_job = select_next_job();
			if ( ! larval_job ) break; // or if we discover there are no jobs to run, then quit

			// ask the job queen to mature the job
			JobOP mature_job = job_queen_->mature_larval_job( larval_job );
			if ( ! mature_job ) {
				// signal from job_queen_ that the inputs for the job are bad
				purge_similar_jobs_which_have_bad_inputs( larval_job );
				continue;
			}

			// run the job
			JobResultOP result;
			try {
				result = mature_job->run();
			} catch ( utility::excn::EXCN_Base const & exception ) {
				process_exception_from_job( larval_job, exception );
			}
			process_job_result( larval_job, result );

		}
		note_round_completed();

	} while ( another_round_remains() );
}

void JobDistributor::set_job_queen( JobQueenOP job_queen )
{
	job_queen_ = job_queen;
}

/// @brief Store the list of jobs to be run in the current round
void JobDistributor::store_jobs_for_current_round( LarvalJobs const & job_list )
{
	jobs_for_current_round_.clear();
	jobs_for_current_round_.resize( job_list.size() );
	std::copy( job_list.begin(), job_list.end(), jobs_for_current_round_.begin() );
}


/// @brief Access to the JobQueen object for derived JobDistributors
JobQueen &
JobDistributor::job_queen() {
	return *job_queen_;
}

LarvalJobs
JobDistributor::determine_jobs_for_next_round() {
	LarvalJobs larval_jobs = job_queen_->determine_job_list();
	njobs_for_round_ = larval_jobs.size();
	next_job_index_  = 1;
	return larval_jobs;
}

bool
JobDistributor::more_jobs_in_current_round() {
	return next_job_index_ <= njobs_for_round_;
}

LarvalJobOP
JobDistributor::select_next_job() {
	while ( next_job_index_ <= njobs_for_round_ ) {
		LarvalJobOP next_job = jobs_for_current_round_[ next_job_index_ ];
		++next_job_index_;
		if ( next_job->bad() ) continue;

		// determine if this job has already been run
		if ( job_queen_->has_job_completed( next_job ) ) {
			continue;
		}
		// ok -- run this job
		job_queen_->mark_job_as_having_begun( next_job );
		return next_job;
	}

	// return 0 if there are no jobs left to run
	return LarvalJobOP( 0 );
}

void
JobDistributor::purge_similar_jobs_which_have_bad_inputs( LarvalJobOP job )
{
	mark_similar_jobs_which_have_bad_inputs_in_job_list( job );
}

void
JobDistributor::mark_similar_jobs_which_have_bad_inputs_in_job_list( LarvalJobOP job )
{
	job_queen_->note_job_completed( job, jd3_job_status_inputs_were_bad );
	for ( core::Size ii = 1; ii <= jobs_for_current_round_.size(); ++ii ) {
		if ( jobs_for_current_round_[ ii ]->completed() || jobs_for_current_round_[ ii ]->bad() ) continue;
		if ( *job == *jobs_for_current_round_[ ii ] ) {
			jobs_for_current_round_[ii]->bad( true );
			job_queen_->note_job_completed( jobs_for_current_round_[ ii ], jd3_job_status_inputs_were_bad );
		}
	}
}

void
JobDistributor::process_exception_from_job( LarvalJobOP job, utility::excn::EXCN_Base const & exception )
{
	TR << "Exception raised by job " << job->job_tag() << " #" << job->nstruct_index() << ":" << std::endl;
	TR << exception.msg() << std::endl;
	job_queen_->note_job_completed( job, jd3_job_status_failed_w_exception );
}


void
JobDistributor::process_job_result( LarvalJobOP job, JobResultOP result ) {
	if ( result->status() == jd3_job_status_success ) {
		job_queen_->completed_job_result( job, result );
	} else if ( result->status() == jd3_job_status_failed_retry ) {
		--next_job_index_;
		// ++count_retries_ ??
	} else if ( result->status() == jd3_job_status_inputs_were_bad ) {
		purge_similar_jobs_which_have_bad_inputs( job );
	} else {
		// failed max retries & failed do not retry
		job_queen_->note_job_completed( job, result->status() );
	}
}


void
JobDistributor::note_round_completed() {}

bool
JobDistributor::another_round_remains()
{
	return job_queen_->more_jobs_remain();
}


}//jd2
}//protocols

