// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/LarvalJob.hh
/// @brief  The definition for clas protocols::jd3::LarvalJob
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_LarvalJob_hh
#define INCLUDED_protocols_jd3_LarvalJob_hh

// Unit headers
#include <protocols/jd3/LarvalJob.fwd.hh>

// Package headers
#include <protocols/jd3/CompletedJobOutput.fwd.hh>
#include <protocols/jd3/InnerLarvalJob.fwd.hh>
#include <protocols/jd3/JobOutputIndex.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <core/types.hh>

// C++ headers
#include <list>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {

/// @brief %protocols::jd3::LarvalJob, which is initialized by the JobQueen, during its
/// intitialize_job_list method, is an immature form of the.  It is matured into a Job
/// that can be run by the JobQueen's mature_job method. A LarvalJob may be one of several
/// LarvalJobs that have the same inputs, but differ in their random number seed (i.e.
/// replicates); two jobs that share the same input are meant to point to the same
/// "InnerLarvalJob" to avoid spending too much memory representing identical data.
///
/// @details Note that in JD2, the currently running Job was globally accessible.  That is not
/// the case in JD3.  Any data that should be output by a running job must be tucked
/// by a Mover inside the JobOutputData to be eventually retrieved by the JobQueen
/// or by a JobDataOutputter at the end of execution.
class LarvalJob : public utility::pointer::ReferenceCount
{
public:


	LarvalJob( InnerLarvalJobOP inner_job, core::Size nstruct_index, core::Size job_index );
	virtual ~LarvalJob() override;

	/// @brief Are these two jobs equivalent? (It does not ask if they
	/// are equal).  Used for the sake of identifying bad jobs. This
	/// method ignores the nstruct_index_, which is unlikely to be equal,
	/// and should instead focus on the data that determines how a job
	/// would run.  If the two jobs share identical input files and parameters
	/// (besides the random number seed), then they should be considered
	/// equivalent.
	/// Derived LarvalJob objects should invoke the base class's version
	/// of this function.
	virtual
	bool
	operator == ( LarvalJob const & rhs ) const;

	virtual
	bool
	operator != ( LarvalJob const & rhs ) const;

	virtual
	void
	show( std::ostream & out ) const;

	/// @brief read access to the inner-job
	InnerLarvalJobCOP inner_job() const;

	/// @brief write access to the inner-job; this should be reserved for the JobQueen only
	/// as she edits the inner job during its construction
	InnerLarvalJobOP nonconst_inner_job() const;

	/// @brief The input tag (a short string, generally), is used to specify the input structure,
	/// but is not a complete description of the LarvalJob, and certainly not the identifier with which
	/// to identify output structures.
	std::string input_tag() const;

	/// @brief The job tag is a combination of the input tag and any other data that the
	/// JobQueen uses to describe the job, or that has been provided by the user in the job-definition
	/// XML file.
	std::string job_tag() const;

	/// @brief For output purposes, construct a string for a particular output of this job given
	/// a JobOutputIndex.
	std::string job_tag_with_index_suffix( JobOutputIndex const & output_index ) const;

	/// @brief The index used to identify which job this is out of many that have identical inputs
	/// but different random number seeds (controlled by the command-line flag "nstruct")
	core::Size nstruct_index() const;

	/// @brief The total number of jobs with the same inputs, but different random number seeds.
	core::Size nstruct_max() const;

	/// @brief The "global" index for this job among all jobs created by the JobQueen
	core::Size job_index() const;

	/// @brief The list of the JobResults required to mature this %LarvalJob, by global index of the
	/// already-executed (Lavral)Jobs and the result_index for that job.
	utility::vector1< JobOutputID > const &
	input_job_result_indices() const;

	void set_status_prefix( std::string const & prefix );
	void set_status_suffix( std::string const & suffix );

	std::string const & status_prefix() const;
	std::string const & status_suffix() const;

	bool completed() const;
	bool to_do() const;
	bool bad() const;
	bool defines_retry_limit() const;
	core::Size retry_limit() const;

	void completed( bool setting );
	void bad( bool setting );
	void retry_limit( core::Size setting );

private:

	/// @brief a pointer to the "heavy" InnerLarvalJob which maintains the data in common with all of the
	/// LarvalJobs (that themselves differ only in their nstruct index).
	InnerLarvalJobOP inner_job_;

	/// @brief which nstruct is this for the given InnerLarvalJob
	core::Size nstruct_index_;

	/// @brief What is the unique global index for this job?  Every %LarvalJob created by the JobQueen
	/// must be assigned a unique index.
	core::Size job_index_;

	/// @brief string giving a brief indication for whether or not the job has failed (but for when the structure
	/// generated by this job should none the less be written as output).
	std::string status_prefix_;

	/// @brief string giving a brief indication for whether or not the job has failed (but for when the structure
	/// generated by this job should none the less be written as output).
	std::string status_suffix_;

	bool completed_;
	bool bad_;
	bool defines_retry_limit_;
	core::Size retry_limit_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	LarvalJob();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // LarvalJob

std::ostream &
operator << ( std::ostream & out, const LarvalJob & job );



} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_LarvalJob )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_LarvalJob_HH
