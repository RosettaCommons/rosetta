// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/LarvalJob.cc
/// @brief  The definition for class protocols::jd3::LarvalJob's methods
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/LarvalJob.hh>

// Package headers
#include <protocols/jd3/InnerLarvalJob.hh>


namespace protocols {
namespace jd3 {

LarvalJob::LarvalJob( InnerLarvalJobOP inner_job, core::Size nstruct_index ) :
	inner_job_( inner_job ),
	nstruct_index_( nstruct_index ),
	completed_( false ),
	bad_( false )
{}

LarvalJob::~LarvalJob() {}

/// @details Are two larval jobs equivalent? This ignores the data that
/// describes the status of the job's execution.  Intentionally ignored
/// data members include:
/// nstruct_index_,
/// status_prefix_,
/// completed_, and
/// bad_
bool
LarvalJob::operator == ( LarvalJob const & rhs ) const {
	return ( inner_job_ == rhs.inner_job_ || ( inner_job_ != 0 && rhs.inner_job_ != 0 && *inner_job_ == *rhs.inner_job_ ));
}

bool
LarvalJob::operator != ( LarvalJob const & rhs ) const {
	return ! ( *this == rhs );
}

void
LarvalJob::show( std::ostream & out ) const {
	out << "LarvalJob show has been stubbed out" << std::endl;
}

/// @brief read access to the inner-job
InnerLarvalJobCOP
LarvalJob::inner_job() const {
	return inner_job_;
}

/// @brief write access to the inner-job; this should be reserved for the JobQueen only
/// as she edits the inner job during its construction
InnerLarvalJobOP
LarvalJob::nonconst_inner_job() const
{
	return inner_job_;
}

/// @brief The input tag (a short string, generally), is used to specify the input structure,
/// but is not a complete description of the LarvalJob, and certainly not the identifier with which
/// to identify output structures.
std::string LarvalJob::input_tag() const {
	return inner_job_->input_tag();
}

/// @brief The job tag is a combination of the input tag and any other data that the
/// JobQueen uses to describe the job, or that has been provided by the user in the job-definition
/// XML file.
std::string LarvalJob::job_tag() const {
	return inner_job_->job_tag();
}

/// @brief The index used to identify which job this is out of many that have identical inputs
/// but different random number seeds (controlled by the command-line flag "nstruct")
core::Size LarvalJob::nstruct_index() const {
	return nstruct_index_;
}

/// @brief The total number of jobs with the same inputs, but different random number seeds.
core::Size LarvalJob::nstruct_max() const {
	return inner_job_->nstruct_max();
}

void LarvalJob::set_status_prefix( std::string prefix ) {
	status_prefix_ = prefix;
}

std::string const &
LarvalJob::status_prefix() const {
	return status_prefix_;
}

bool LarvalJob::completed() const {
	return completed_;
}

bool LarvalJob::to_do() const {
	return !completed_ && !bad_;
}

bool LarvalJob::bad() const {
	return bad_;
}

void LarvalJob::completed( bool setting ) {
	completed_ = setting;
}

void LarvalJob::bad( bool setting ) {
	bad_ = setting;
}

std::ostream &
operator << ( std::ostream & out, const LarvalJob & job )
{
	job.show( out );
	return out;
}

} // namespace jd3
} // namespace protocols
