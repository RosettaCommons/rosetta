// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/LarvalJob.cc
/// @brief  The definition for class protocols::jd3::LarvalJob's methods
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/LarvalJob.hh>

// Package headers
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/JobOutputIndex.hh>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <cmath>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {

LarvalJob::LarvalJob( InnerLarvalJobOP inner_job, core::Size nstruct_index, core::Size job_index ) :
	inner_job_(std::move( inner_job )),
	nstruct_index_( nstruct_index ),
	job_index_( job_index ),
	completed_( false ),
	bad_( false ),
	defines_retry_limit_( false ),
	retry_limit_( 1 )
{}

LarvalJob::~LarvalJob() = default;

/// @details Are two larval jobs equivalent? This ignores the data that
/// describes the status of the job's execution.  Intentionally ignored
/// data members include:
/// nstruct_index_,
/// status_prefix_,
/// completed_, and
/// bad_
bool
LarvalJob::operator == ( LarvalJob const & rhs ) const {
	return ( inner_job_ == rhs.inner_job_ || ( inner_job_ != nullptr && rhs.inner_job_ != nullptr && *inner_job_ == *rhs.inner_job_ ));
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

/// @details The number of leading zeros is determined by the maximum number of jobs for the inner-job
/// that this job points at.  For 9999  nstruct, there should be 4 digits: 1 + int(log10( 9999 )) = 1 + int( 3.9999 ) = 4
/// For 10K nstruct, there should be 5 digits; 1 + int( log10( 10K )) = 5

std::string LarvalJob::job_tag_with_index_suffix( JobOutputIndex const & output_index ) const
{
	return inner_job_->job_tag() + "_" +
		ObjexxFCL::lead_zero_string_of( output_index.primary_output_index, std::max( 4, 1 + int( std::log10( output_index.n_primary_outputs ))) ) +
		( output_index.secondary_output_index == 1 && output_index.n_secondary_outputs == 1 ?
		"" : "_" + ObjexxFCL::lead_zero_string_of( output_index.secondary_output_index, std::max( 4, 1 + int( std::log10( output_index.n_secondary_outputs ))))) ;


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

core::Size LarvalJob::job_index() const {
	return job_index_;
}

utility::vector1< JobOutputID > const &
LarvalJob::input_job_result_indices() const
{
	return inner_job_->input_job_result_indices();
}

void LarvalJob::set_status_prefix( std::string const & prefix ) {
	status_prefix_ = prefix;
}

void LarvalJob::set_status_suffix( std::string const & suffix )
{
	status_suffix_ = suffix;
}

std::string const &
LarvalJob::status_prefix() const {
	return status_prefix_;
}

std::string const &
LarvalJob::status_suffix() const {
	return status_suffix_;
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

bool LarvalJob::defines_retry_limit() const {
	return defines_retry_limit_;
}

core::Size LarvalJob::retry_limit() const {
	return retry_limit_;
}

void LarvalJob::completed( bool setting ) {
	completed_ = setting;
}

void LarvalJob::bad( bool setting ) {
	bad_ = setting;
}

void LarvalJob::retry_limit( core::Size setting ) {
	defines_retry_limit_ = true;
	retry_limit_ = setting;
}


std::ostream &
operator << ( std::ostream & out, const LarvalJob & job )
{
	job.show( out );
	return out;
}

} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
protocols::jd3::LarvalJob::LarvalJob() : nstruct_index_( 0 ), job_index_( 0 ) {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::LarvalJob::save( Archive & arc ) const {
	arc( CEREAL_NVP( inner_job_ ) ); // InnerLarvalJobOP
	arc( CEREAL_NVP( nstruct_index_ ) ); // core::Size
	arc( CEREAL_NVP( job_index_ ) ); // core::Size
	arc( CEREAL_NVP( status_prefix_ ) ); // std::string
	arc( CEREAL_NVP( status_suffix_ ) ); // std::string
	arc( CEREAL_NVP( completed_ ) ); // _Bool
	arc( CEREAL_NVP( bad_ ) ); // _Bool
	arc( CEREAL_NVP( defines_retry_limit_ ) ); // _Bool
	arc( CEREAL_NVP( retry_limit_ ) ); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::LarvalJob::load( Archive & arc ) {
	arc( inner_job_ ); // InnerLarvalJobOP
	arc( nstruct_index_ ); // core::Size
	arc( job_index_ ); // core::Size
	arc( status_prefix_ ); // std::string
	arc( status_suffix_ ); // std::string
	arc( completed_ ); // _Bool
	arc( bad_ ); // _Bool
	arc( defines_retry_limit_ ); // _Bool
	arc( retry_limit_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::LarvalJob );
CEREAL_REGISTER_TYPE( protocols::jd3::LarvalJob )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_LarvalJob )
#endif // SERIALIZATION
