// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/JobsContainer.cc
/// @brief  A container for an array of JobOPs.
/// @details By making this its own class, it's easier to update the list if, for example, we
///      don't want to store the whole list in memory, but generate it on-the-fly,
///      piecewise.
/// @author  Vikram K. Mulligan, Baker Laboratory (vmullig@uw.edu)

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

// Unit headers
#include <protocols/jd2/JobsContainer.hh>

// Project headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/Job.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/sys_util.hh>
#include <numeric/random/random_permutation.hh>

// C++ headers
#include <string>
#include <ctime>

#include <utility/vector1.hh>


static basic::Tracer TR( "protocols.jd2.JobsContainer" );

namespace protocols {
namespace jd2 {

JobsContainer::JobsContainer() :
	joblist_(),
	highest_job_index_(0),
	total_jobs_(0),
	total_jobs_set_(false),
	job_inputter_(), //NULL initially
	force_job_purging_(false)
{
}

JobsContainer::~JobsContainer()= default;

/// @brief Get a specific job, by number.
/// @details Should work even if jobs have been deleted, since this uses a
/// map instead of an array.
JobOP JobsContainer::operator [](core::Size const index) {

	if ( index < 1 || index>total_jobs_ ) {
		utility_exit_with_message(
			"In protocols::jd2::JobsContainer [] operator: index is out of range (less than 1 or greater than total jobs).\n"
		);
	}

	//if(index <= highest_job_index() && !has_job(index) && TR.visible()) { TR << "Warning!  Already-deleted job " << index << " requested!" << std::endl; TR.flush(); } //DELETE ME

	if ( index>highest_job_index() ) {
		if ( !job_inputter_ ) {
			utility_exit_with_message(
				"In protocols::jd2::JobsContainer [] operator: NULL owning pointer for the job inputter provided."
			);
		}
		job_inputter_->update_jobs_list( get_self_ptr() ); //Purge old jobs and add more.
		if ( index>highest_job_index() ) { //If we're still asking for a job greater than the highest in memory...
			bool success(false);
			if ( force_job_purging() ) {
				if ( TR.visible() ) {
					TR << "Forcing job purging..." << std::endl;
					TR.flush();
				}
				core::Size purge_iter(0);
				while ( highest_job_index() < size() && success==false ) { //While we have additional jobs that we could load into memory
					if ( TR.visible() ) {
						TR << "Purge iteration " << ++purge_iter << std::endl;
						TR.flush();
					}
					for ( core::Size i=1, imax=highest_job_index(); i<=imax; ++i ) { //Set everything as deletable
						if ( has_job(i) ) joblist_[i]->set_can_be_deleted(true);
					}
					job_inputter_->update_jobs_list( get_self_ptr() ); //Repopulate the list.
					if ( index <= highest_job_index() ) { //Is the index now within the list bounds?
						success=true;
						break;
					}
				}
			}
			if ( !success ) {
#ifdef USEMPI
							utility_exit_with_message(
								"In protocols::jd2::JobsContainer [] operator: could not update the jobs list to get requested job.  This could be because the maximum number of jobs was already present in memory.  (Note that this can happen in MPI mode if the number of processes is set higher than the maximum number of jobs that can be held in memory at any given time.  If this is an MPI job, please check that the value specified with the -jd2:max_nstruct_in_memory option is greater than the number of processes specified with mpirun.)\n"
							);
#else
				utility_exit_with_message(
					"In protocols::jd2::JobsContainer [] operator: could not update the jobs list to get requested job.  This could be because the maximum number of jobs was already present in memory.\n"
				);
#endif
			}
		}
	}

	return joblist_.at(index);
}

/// @brief Get a specific job, by number (const-access).
/// @details Should work even if jobs have been deleted, since this uses a
/// map instead of an array.
JobCOP JobsContainer::operator [](core::Size const index) const {
	if ( index < 1 || index>total_jobs_ ) {
		utility_exit_with_message(
			"In protocols::jd2::JobsContainer [] const operator: index is out of range (less than 1 or greater than total jobs).\n"
		);
	}
	if ( index>highest_job_index() ) {
		utility_exit_with_message(
			"In protocols::jd2::JobsContainer [] const operator: index is out of range of jobs currently loaded in memory.  Since this is const-access to the jobs list, we cannot load more jobs.\n"
		);
	}
	return joblist_.at(index);
}

/// @brief Assignment operator.
///
JobsContainer &
JobsContainer::operator=( JobsContainer const & src ) {
	joblist_ = src.joblist_;
	highest_job_index_=src.highest_job_index_;
	total_jobs_=src.total_jobs_;
	total_jobs_set_=src.total_jobs_set_;
	job_inputter_=src.job_inputter_; //Do not clone!  Use OP directly!
	force_job_purging_=src.force_job_purging_;
	return *this;
}

/// @brief Add a job to the list of jobs.
///
void JobsContainer::push_back( JobOP new_job ) {
	++highest_job_index_;
	joblist_.insert( std::pair<core::Size, JobOP>(highest_job_index_, new_job) );
	if ( !total_jobs_set_ ) total_jobs_=highest_job_index_;
	else {
		if ( highest_job_index_ > total_jobs_ ) {
			utility_exit_with_message( "In protocols::jd2::JobsContainer::push_back(): The number of jobs was set, but more jobs were added to the container than allowed.\n" );
		}
	}
	return;
}

/// @brief Set the total number of jobs.
/// @details This overrides whatever the length of the joblist_ is.
void JobsContainer::set_total_jobs( core::Size const total_jobs ) {
	if ( highest_job_index_ > total_jobs ) {
		utility_exit_with_message( "In protocols::jd2::JobsContainer::set_total_jobs(): The number of jobs was set to a lower number than the number already allocated.\n" );
	}
	total_jobs_=total_jobs;
	total_jobs_set_=true;
	return;
}

/// @brief Clear the jobs list
///
void JobsContainer::clear() {
	joblist_.clear();
	highest_job_index_=0;
	total_jobs_=0;
	total_jobs_set_=false;
	return;
}

/// @brief Access the last element.
///
JobOP JobsContainer::back() {
	return joblist_[highest_job_index_];
}

/// @brief Randomize the order of elements (the map keys)
///
void JobsContainer::shuffle() {
	utility::vector1<core::Size> keys;
	for ( auto & it : joblist_ ) {
		keys.push_back( it.first );
	}
	numeric::random::random_permutation( keys, numeric::random::rg() );
	std::map<core::Size, JobOP> newjoblist_;
	core::Size j=1;
	for ( auto & it : joblist_ ) {
		newjoblist_.insert( std::pair<core::Size, JobOP>( keys[j], it.second ) );
		++j;
	}
	joblist_=newjoblist_;
	return;
}

/// @brief Erase an element in the jobs list
void JobsContainer::erase( core::Size const index ) {
	//if(index==highest_job_index_) --highest_job_index_;
	std::map< core::Size, JobOP >::iterator it;
	it=joblist_.find(index);
	joblist_.erase( it );
	return;
}

/// @brief Does the job with the given index exist in the currently-loaded list of jobs?
///
bool JobsContainer::has_job( core::Size const index ) const {
	return (joblist_.count( index ) > 0);
}

/// @brief Can the job with the given index be deleted?
///
bool JobsContainer::can_be_deleted ( core::Size const index ) const {
	return ( joblist_.at( index )->can_be_deleted() );
}

/// @brief Get a list of job indices currently in memory.
/// @details The output vector is cleared and populated with the current job indices in memory.
void JobsContainer::get_loaded_job_indices( utility::vector1 < core::Size > &output) const {
	output.clear();
	output.reserve( joblist_.size() );
	for ( auto const & it : joblist_ ) {
		output.push_back( it.first );
	}
	return;
}

/// @brief Mark all jobs currently in memory as deletable.
/// @details This will result in all jobs being purged the next time a higher-index job is requeted.
void JobsContainer::set_all_current_jobs_as_deletable() {

	return;
}

} // jd2
} // protocols
