// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/JobTracker.hh
/// @brief A simple class for tracking job progress within JD3.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_jd3_JobTracker_hh
#define INCLUDED_protocols_jd3_JobTracker_hh

#include <protocols/jd3/JobTracker.fwd.hh>
#include <protocols/jd3/JobQueen.fwd.hh>

// Core headers
#include <core/types.hh>

// Protocol headers
#include <protocols/jd3/LarvalJob.fwd.hh>
#include <protocols/jd3/Job.fwd.hh>

// numeric headers
#include <numeric/DiscreteIntervalEncodingTree.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <unordered_map>
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {

/// @brief A simple class for tracking job progress within JD3.
class JobTracker : public utility::pointer::ReferenceCount {

public:

	JobTracker();

	JobTrackerOP
	clone() const;

public:

	///@brief Note the completed job by Job Dag node and by Status.
	/// Only the JobQueen is allowed to call this method.
	void
	track_completed_job( JQKey key, LarvalJob const & larval_job, JobStatus status);

	void
	track_starting_job_list( JQKey key, LarvalJobs const & starting_jobs);

	///@brief Increment the current job index tracked by the JobTracker.
	///
	///@details Used to set the job_index for LarvalJobs during determine_job_list
	///
	void
	increment_current_job_index();

public:

	/// @brief Read access for jobs have been given out to the JD through determine_job_list.
	numeric::DiscreteIntervalEncodingTree< core::Size > const &
	started_jobs() const;

	/// @brief Read access to all job indexes started for a particular job node.
	numeric::DiscreteIntervalEncodingTree< core::Size > const &
	started_jobs_for_node( core::Size job_dag_node ) const;

	///@brief Last job index of a particular pose id from all starting jobs.
	core::Size
	last_job_for_input_source_id( core::Size input_pose_index ) const;

	///@brief Get the map of last job indexes and input source ids
	std::map< core::Size, core::Size > const &
	last_jobs_for_inputs_sources() const;

public:
	/// @brief Read access for which jobs have completed and how; if a job-id is a member
	/// of this DIET, then it has completed (either in success or failure).
	numeric::DiscreteIntervalEncodingTree< core::Size > const &
	completed_jobs() const;

	numeric::DiscreteIntervalEncodingTree< core::Size > const &
	completed_jobs_for_node( core::Size job_dag_node ) const;

	/// @brief Read access for which jobs have completed and how; if a job-id is a member
	/// of this DIET, then it completed successfully.
	numeric::DiscreteIntervalEncodingTree< core::Size > const &
	successful_jobs() const;

	/// @brief Read access for which jobs have completed and how; if a job-id is a member
	/// of this DIET, then it completed unsuccessfully.
	numeric::DiscreteIntervalEncodingTree< core::Size > const &
	failed_jobs() const;

	/// @brief Read access for which jobs have completed and how; if a job-id is a member
	/// of this DIET, then the job has had all of its results output or discarded.
	//numeric::DiscreteIntervalEncodingTree< core::Size > const &
	//processed_jobs() const;

public:

	///@brief Get the current larval job index
	core::Size
	current_job_index() const;

private:
	typedef numeric::DiscreteIntervalEncodingTree< core::Size > SizeDIET;

	std::unordered_map< core::Size, SizeDIET> completed_jobs_by_dag_node_;
	std::unordered_map< core::Size, SizeDIET> started_jobs_by_dag_node_;

	SizeDIET completed_jobs_;
	SizeDIET started_jobs_;

	SizeDIET successful_jobs_;
	SizeDIET failed_jobs_;
	SizeDIET previously_completed_jobs_;

	std::map< core::Size, core::Size > last_job_for_input_source_;

	// The index of the last larval job that we have created. Incremented within
	// expand_job_list().
	Size current_global_job_index_ = 0;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //protocols
} //jd3



#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_JobTracker )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_JobTracker_hh





