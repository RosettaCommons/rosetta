// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/JobTracker.cc
/// @brief A simple class for tracking job progress within JD3.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/jd3/JobTracker.hh>

// Protocol headers
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/InputSource.hh>

// Numeric headers
#include <numeric/DiscreteIntervalEncodingTree.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>



static basic::Tracer TR( "protocols.jd3.JobTracker" );


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
//#include <cereal/types/hash.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {

JobTracker::JobTracker():
	utility::pointer::ReferenceCount()
{
	completed_jobs_by_dag_node_.max_load_factor( 0.7 );
	started_jobs_by_dag_node_.max_load_factor( 0.7 );
}


JobTrackerOP
JobTracker::clone() const {
	return utility::pointer::make_shared< JobTracker>( *this );
}



void
JobTracker::track_starting_job_list( JQKey , LarvalJobs const & larval_jobs){
	for ( auto const & larval_job : larval_jobs ) {
		core::Size job_dag_node = larval_job->job_node();
		core::Size job_id = larval_job->job_index();

		started_jobs_.insert( job_id );
		if ( ! started_jobs_by_dag_node_.count(job_dag_node) ) {
			numeric::DiscreteIntervalEncodingTree< core::Size > intervals;
			intervals.insert(job_id);
			started_jobs_by_dag_node_[job_dag_node] = intervals;
		} else {
			started_jobs_by_dag_node_[job_dag_node].insert( job_id );
		}

		//Only track if we actually have input sources to track.
		if ( larval_job->inner_job()->n_input_sources() > 0 ) {
			for ( core::Size source = 1; source <= larval_job->inner_job()->n_input_sources(); ++source ) {
				core::Size source_id = larval_job->inner_job()->input_source(source).source_id();

				//If we don't have a PoseID set (ie 0), don't track it.
				if ( (source_id != 0) && (!last_job_for_input_source_.count( source_id )) ) {
					last_job_for_input_source_[ source_id ] = larval_job->job_index();
				} else if ( (source_id !=0)  && (last_job_for_input_source_[ source_id ] < larval_job->job_index()) ) {
					last_job_for_input_source_[ source_id ] = larval_job->job_index();
				}
			}
		}
	}
}

/// @brief Read access for jobs have been given out to the JD through determine_job_list.
numeric::DiscreteIntervalEncodingTree< core::Size > const &
JobTracker::started_jobs() const{
	return started_jobs_;
}

/// @brief Read access to all job indexes started for a particular job node.
numeric::DiscreteIntervalEncodingTree< core::Size > const &
JobTracker::started_jobs_for_node( core::Size job_dag_node ) const{
	return started_jobs_by_dag_node_.at(job_dag_node);
}

///@brief Last job index of a particular pose id from all starting jobs.
core::Size
JobTracker::last_job_for_input_source_id( core::Size input_source_index ) const{
	return last_job_for_input_source_.at( input_source_index);
}

std::map< core::Size, core::Size > const &
JobTracker::last_jobs_for_inputs_sources() const {
	return last_job_for_input_source_;
}

///@brief Note the completed job by Job Dag node and by Status.
void
JobTracker::track_completed_job( JQKey, LarvalJob const & larval_job, JobStatus status){
	core::Size job_dag_node = larval_job.job_node();
	core::Size job_id = larval_job.job_index();

	completed_jobs_.insert( job_id );
	if ( ! completed_jobs_by_dag_node_.count(job_dag_node) ) {
		numeric::DiscreteIntervalEncodingTree< core::Size > intervals;
		intervals.insert(job_id);
		completed_jobs_by_dag_node_[job_dag_node] = intervals;
	} else {
		completed_jobs_by_dag_node_[job_dag_node].insert(job_id);
	}

	if ( status == jd3_job_status_success ) {
		successful_jobs_.insert( job_id );
	} else if ( status == jd3_job_previously_executed ) {
		previously_completed_jobs_.insert( job_id );
	} else {
		failed_jobs_.insert( job_id );
	}

	//for ( auto it = last_job_for_input_source_.begin(); it != last_job_for_input_source_.end(); /*no inc*/ ) {
	//
	// if ( it->second == job_id ) {
	//  it = last_job_for_input_pose_.erase( it );
	// } else {
	//  ++it;
	// }
	//}

}


/// @brief Read access for which jobs have completed and how; if a job-id is a member
/// of this DIET, then it has completed (either in success or failure).
numeric::DiscreteIntervalEncodingTree< core::Size > const &
JobTracker::completed_jobs() const{
	return completed_jobs_;
}

numeric::DiscreteIntervalEncodingTree< core::Size > const &
JobTracker::completed_jobs_for_node( core::Size job_dag_node ) const{
	return completed_jobs_by_dag_node_.at(job_dag_node);
}

/// @brief Read access for which jobs have completed and how; if a job-id is a member
/// of this DIET, then it completed successfully.
numeric::DiscreteIntervalEncodingTree< core::Size > const &
JobTracker::successful_jobs() const{
	return successful_jobs_;
}

/// @brief Read access for which jobs have completed and how; if a job-id is a member
/// of this DIET, then it completed unsuccessfully.
numeric::DiscreteIntervalEncodingTree< core::Size > const &
JobTracker::failed_jobs() const{
	return failed_jobs_;
}

void
JobTracker::increment_current_job_index(){
	current_global_job_index_+=1;
}

///@brief Get read access to the current job index
core::Size
JobTracker::current_job_index() const {
	return current_global_job_index_;
}



} //protocols
} //jd3







#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::JobTracker::save( Archive & arc ) const {
	arc( CEREAL_NVP( completed_jobs_by_dag_node_ ) ); // std::unordered_map<core::Size, SizeDIET>
	arc( CEREAL_NVP( started_jobs_by_dag_node_ ) ); // std::unordered_map<core::Size, SizeDIET>
	arc( CEREAL_NVP( completed_jobs_ ) ); // SizeDIET
	arc( CEREAL_NVP( started_jobs_ ) ); // SizeDIET
	arc( CEREAL_NVP( successful_jobs_ ) ); // SizeDIET
	arc( CEREAL_NVP( failed_jobs_ ) ); // SizeDIET
	arc( CEREAL_NVP( previously_completed_jobs_ ) ); // SizeDIET
	arc( CEREAL_NVP( last_job_for_input_source_ ) ); // std::map<core::Size, core::Size>
	arc( CEREAL_NVP( current_global_job_index_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::JobTracker::load( Archive & arc ) {
	arc( completed_jobs_by_dag_node_ ); // std::unordered_map<core::Size, SizeDIET>
	arc( started_jobs_by_dag_node_ ); // std::unordered_map<core::Size, SizeDIET>
	arc( completed_jobs_ ); // SizeDIET
	arc( started_jobs_ ); // SizeDIET
	arc( successful_jobs_ ); // SizeDIET
	arc( failed_jobs_ ); // SizeDIET
	arc( previously_completed_jobs_ ); // SizeDIET
	arc( last_job_for_input_source_ ); // std::map<core::Size, core::Size>
	arc( current_global_job_index_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::JobTracker );
CEREAL_REGISTER_TYPE( protocols::jd3::JobTracker )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_JobTracker )
#endif // SERIALIZATION
