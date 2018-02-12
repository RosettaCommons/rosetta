// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/dag_node_managers/NodeManager.hh
/// @brief Base class for the family of JD3 Node Managers. This class is intended to offload some of the result-sorting logic from a Job Queen for a single job-dag node.

/// @detailed The base class is somewhat messy, so there are a few derived classes in protocols/jd3/DerivedNodeManagers.hh that simplify some of these interfaces by specializing for certain cases. Features include:
/// - Sorting results by some metric determined by the user (more negative values are considered "better")
/// - Partitioning results into separate bins
/// - Keeping track of jobs that have been susbmitted and the global job offset
/// - Finishing early (if result_threshold argument in the constructor is != 0) if enough results come in
/// - Determining which jobs results should be discarded and which should be kept
/// - See here for more information: https://www.rosettacommons.org/docs/latest/development_documentation/tutorials/jd3_derived_jq/classes/node_manager
/*!
Underlying data structure resembles a 2D vector.
Each row (outer dimension) represents a different partition and each partition
holds a sorted vector of results that is pruned to make sure it has no more than
the maximum number of elements allowed. elements that are pruned out are added to a list of job results to discard.

To get the most out of this class, call:

1) get_next_local_jobid() during JobQueen::determine_job_list()

2) note_job_completed() during JobQueen::note_job_completed()

3) register_result() during JobQueen::completed_job_summary()

4) append_job_results_that_should_be_discarded() during JobQueen::job_results_that_should_be_discarded()

5) get_nth_job_result_id() during JobQueen::determine_job_list() for downstream job dag nodes
*/
/// @author Jack Maguire, jackmaguire1444@gmail.com


#ifndef INCLUDED_protocols_jd3_dag_node_managers_NodeManager_HH
#define INCLUDED_protocols_jd3_dag_node_managers_NodeManager_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/jd3/dag_node_managers/NodeManager.fwd.hh>
#include <protocols/jd3/dag_node_managers/NodeManagerStorageMatrix.hh>

// Package headers
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/JobResult.hh>
#include <protocols/jd3/JobSummary.hh>


namespace protocols {
namespace jd3 {
namespace dag_node_managers {


class NodeManager : public utility::pointer::ReferenceCount {

public:

	//constructor
	///@param job_offset The node manager can only represent nodes where the jobids form a continuous range. That range should start at job_offset+1.
	///@param num_jobs_total The range mentioned previous should end with job_offset+num_jobs_total (inclusive)
	///@param num_partitions The total number of partitions you want
	///@param num_results_to_keep_for_part maximum number of results you want to keep per partition. num_results_to_keep_for_part.size() should equal num_partitions
	///@param result_threshold_per_part If you want to have a result threshold, define here the threshold for each partition. result_threshold_per_part.size() should equal num_partitions
	NodeManager(
		core::Size job_offset,
		core::Size num_jobs_total,
		core::Size num_partitions,
		utility::vector1< core::Size > num_results_to_keep_for_part,
		utility::vector1< core::Size > result_threshold_per_part = utility::vector1< core::Size > ( 0 ),
		bool return_results_depth_first = false
	);

	//destructor
	~NodeManager() override;

public:

	///@brief please call this from your job queen's note_job_completed() function
	void note_job_completed(
		core::Size global_job_id,
		core::Size nresults
	);

	///@brief insert this result into the sorted container.
	///Please specify which partition you want to put it in if you are not using the SimpleNodeManager
	virtual void register_result(
		core::Size global_job_id,
		core::Size local_result_id,
		core::Real score,
		core::Size partition = 1
	);

	///@brief This class can be used to determine which job should be submitted next.
	/// A value of 0 means that we are done submitting for this dag node
	inline core::Size get_next_local_jobid() {
		if ( done_submitting() ) return 0;
		return ++num_jobs_submitted_;
	}

	///@brief This does not erase the list! We just add job result ids that have been eliminated from the results to keep
	inline void append_job_results_that_should_be_discarded(
		std::list< jd3::JobResultID > & list
	){
		list.splice( list.end(), job_results_that_should_be_discarded_ );
	}

	///@brief this has better future-proofing over results_to_keep()[ int ]. Returns {0,0} if no result
	jd3::JobResultID get_nth_job_result_id( core::Size n );

public://status checkers

	inline bool done_submitting() const {
		return stopped_early_ || num_jobs_submitted_ == num_jobs_total_;
	}

	inline bool jobs_are_still_running(){
		return num_jobs_completed_ < num_jobs_submitted_;
	}

	inline bool all_results_are_in() const {
		return done_submitting()
			&& num_jobs_submitted_ == num_jobs_completed_
			&& num_results_received_ == num_results_total_;
	}

	inline bool all_waste_is_discarded() const {
		return all_results_are_in() && job_results_that_should_be_discarded_.empty();
	}

public: //getters, setters

	///@brief old way to access the results. get_nth_job_result_id() is prefered.
	utility::vector1< result_elements > results_to_keep();

	inline core::Size num_jobs() const {
		return num_jobs_total_;
	}

	inline core::Size num_jobs_submitted() const {
		return num_jobs_submitted_;
	}

	inline core::Size num_jobs_completed() const {
		return num_jobs_completed_;
	}

	inline core::Size num_results_to_keep() const {
		return num_results_to_keep_;
	}

	inline core::Size job_offset() const {
		return job_offset_;
	}

	inline core::Size num_results_total() const {
		return num_results_total_;
	}

	inline void stop_early(){ stopped_early_ = true; }

	inline bool stopped_early() const { return stopped_early_; };

protected:
	inline core::Size num_results_received(){
		return num_results_received_;
	}

	inline bool ready_to_finish_early() const {
		for ( core::Size ii = 1; ii <= num_partitions_; ++ii ) {
			if ( num_results_received_for_part_[ ii ] < result_threshold_per_part_[ ii ] ) {
				return false;
			}
		}
		return true;
	}

	void set_return_results_depth_first( bool setting ){
		results_to_keep_.set_return_results_depth_first( setting );
	}

private:

	core::Size job_offset_;//The first job for this node will be job_offset_ + 1
	core::Size num_jobs_total_;

	core::Size num_results_to_keep_;
	utility::vector1< core::Size > num_results_to_keep_for_part_;

	core::Size num_results_total_;
	core::Size num_results_received_;
	utility::vector1< core::Size > num_results_received_for_part_;

	bool result_threshold_;
	utility::vector1< core::Size > result_threshold_per_part_;
	bool stopped_early_;

	core::Size num_jobs_submitted_;
	core::Size num_jobs_completed_;

	core::Size num_partitions_;
	//utility::vector1< utility::vector1< result_elements > > results_to_keep_;
	NodeManagerStorageMatrix results_to_keep_;

	//mutable utility::vector1< result_elements > all_results_to_keep_;

	std::list< jd3::JobResultID > job_results_that_should_be_discarded_;
};

inline
void
NodeManager::note_job_completed( core::Size, core::Size nresults ){//first arg is global_job_id
	++num_jobs_completed_;
	num_results_total_ += nresults;
}

inline
utility::vector1< result_elements >
NodeManager::results_to_keep() {
	//Have to stop early because this freezes the result matrix.
	//This is one of the main reasons not to use this method.
	stopped_early_ = true;
	return results_to_keep_.linear_vector_of_results();
}


inline
jd3::JobResultID
NodeManager::get_nth_job_result_id( core::Size n ) {
	return results_to_keep_.get_nth_element( n ).job_result_id();
}


} //dag_node_managers
} //jd3
} //protocols

#endif
