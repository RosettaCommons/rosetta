// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/standard/PreliminaryLarvalJobTracker.cc
/// @brief A class that tracks PreliminaryLarvalJobs in the SJQ.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/jd3/standard/PreliminaryLarvalJobTracker.hh>
#include <protocols/jd3/standard/PreliminaryLarvalJob.hh>

//Protocol includes
#include <protocols/jd3/InputSource.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>


#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.jd3.standard.PreliminaryLarvalJobTracker" );


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace standard {

PreliminaryLarvalJobTracker::PreliminaryLarvalJobTracker():
	utility::pointer::ReferenceCount()
{

}


void
PreliminaryLarvalJobTracker::initialize_tracker( utility::vector1< PreliminaryLarvalJob > const & prelim_jobs ){

	if ( prelim_jobs.size() == 0 ) return;

	//debug_assert(prelim_jobs.size() != 0 );

	preliminary_job_node_inds_.clear();
	preliminary_job_nodes_assigned_.clear();
	pjn_job_ind_end_.clear();
	pjn_job_ind_begin_.clear();
	incomplete_pjns_using_input_pose_.clear();
	outstanding_job_count_.clear();
	pjn_max_nstruct_.clear();

	core::Size total_preliminary_job_nodes = prelim_jobs.size();
	preliminary_job_node_inds_.resize( total_preliminary_job_nodes, 0 );
	preliminary_job_nodes_assigned_.resize( total_preliminary_job_nodes, false );
	pjn_job_ind_begin_.resize( total_preliminary_job_nodes, 0 );
	pjn_job_ind_end_.resize( total_preliminary_job_nodes, 0 );

	for ( core::Size ii = 1; ii <= total_preliminary_job_nodes; ++ii ) {
		preliminary_job_node_inds_[ ii ] = ii;
		pjn_job_ind_begin_[ ii ] = pjn_job_ind_end_[ ii ] = 0;
		preliminary_job_nodes_assigned_[ ii ] = 0;
	}

	//Initialize per-preliminary job tracking data.
	for ( core::Size job_node = 1; job_node <= prelim_jobs.size(); ++job_node ) {

		PreliminaryLarvalJob const & prelim_job = prelim_jobs[job_node];

		core::Size nstruct = prelim_job.inner_job->nstruct_max();
		core::Size pose_id = prelim_job.inner_job->input_source().source_id();
		outstanding_job_count_.push_back( nstruct );
		pjn_max_nstruct_.push_back( nstruct );
		incomplete_pjns_using_input_pose_[ pose_id ].insert(job_node);

		if ( ! non_deallocated_input_poses_.contains(pose_id) ) {
			non_deallocated_input_poses_.push_back(pose_id);
		}
	}
	TR.Debug << "Preliminary Job Nodes of size " << total_preliminary_job_nodes << " initialized" << std::endl;
}

void
PreliminaryLarvalJobTracker::track_job_completed(LarvalJobCOP job ){
	core::Size pjn_index = get_job_node_for_job_index( job->job_index() );
	if ( pjn_index != 0 ) {
		--outstanding_job_count_[ pjn_index ];
		TR.Debug << "Outstanding job_count for PJN "<< pjn_index << " " << outstanding_job_count_[ pjn_index ] << std::endl;

		// Now we need to note that this PJN is no longer blocking the deallocation of
		// any of the input poses that we might still be holding on to
		if ( get_job_node_complete(pjn_index) ) {
			TR.Debug << "Job Node Complete: " << pjn_index << std::endl;
			for ( core::Size ii = 1; ii <= job->inner_job()->n_input_sources(); ++ii ) {
				InputSource const & ii_source = job->inner_job()->input_source( ii );
				incomplete_pjns_using_input_pose_[ ii_source.source_id() ].erase( pjn_index );
				TR.Debug << "Deallocating source:"  << ii_source.source_id() << " pjn: " << pjn_index << std::endl;
			}

		}
	}
}

void
PreliminaryLarvalJobTracker::track_job_node_assigned( core::Size job_node_index, core::Size last_global_nstruct_index_for_node){
	pjn_job_ind_end_[ job_node_index] = last_global_nstruct_index_for_node;
	preliminary_job_nodes_assigned_[ job_node_index] = true;
}

void
PreliminaryLarvalJobTracker::track_job_node_being_assigned(core::Size job_node_index, core::Size global_start_index, core::Size global_end_index){

	if ( pjn_job_ind_begin_[ job_node_index ] == 0 ) {
		pjn_job_ind_begin_[ job_node_index ] = global_start_index;
	}
	pjn_job_ind_end_[ job_node_index ] = global_end_index;
}

bool
PreliminaryLarvalJobTracker::get_job_node_complete( core::Size job_node_index ) const{
	return outstanding_job_count_[ job_node_index ] == 0;
}

bool
PreliminaryLarvalJobTracker::get_job_node_assigned( core::Size job_node_index ) const {
	return preliminary_job_nodes_assigned_[ job_node_index ];
}

core::Size
PreliminaryLarvalJobTracker::get_job_node_for_job_index( core::Size job_id ) const {
	// TO DO: Replace with binary search?
	for ( core::Size ii = 1; ii <= pjn_job_ind_begin_.size(); ++ii ) {
		if ( job_id >= pjn_job_ind_begin_[ ii ] && job_id <= pjn_job_ind_end_[ ii ] ) {
			return ii;
		}
	}
	return 0;
}

bool
PreliminaryLarvalJobTracker::input_pose_no_longer_needed(core::Size pose_index) const {
	return ( incomplete_pjns_using_input_pose_.at( pose_index ).empty() );
}

///@brief Get the starting global job index that starts for a particular preliminary job node
core::Size
PreliminaryLarvalJobTracker::get_job_index_starting_job_node(core::Size job_node){
	if ( job_node > pjn_job_ind_begin_.size() ) {
		return 0;
	} else {
		return pjn_job_ind_begin_[job_node];
	}
}

///@brief Get the ending global job index that starts for a particular preliminary job node
core::Size
PreliminaryLarvalJobTracker::get_job_index_ending_job_node(core::Size job_node){
	if ( job_node > pjn_job_ind_end_.size() ) {
		return 0;
	} else {
		return pjn_job_ind_end_[job_node];
	}
}

utility::vector1< core::Size >
PreliminaryLarvalJobTracker::get_preliminary_job_node_indices() const {
	return preliminary_job_node_inds_;
}

///@brief Track a deallocated input pose. (indicate that this pose id has been deallocated)
void
PreliminaryLarvalJobTracker::track_deallocated_input_pose( core::Size pose_id ){
	non_deallocated_input_poses_ = non_deallocated_input_poses_.pop(pose_id);
	incomplete_pjns_using_input_pose_[pose_id].clear(); //Make sure this is cleared if for some reason we need to deallocate it.
}

utility::vector1< core::Size >
PreliminaryLarvalJobTracker::get_input_poses_to_deallocate() const{

	utility::vector1< core::Size > final_deallocation_ids;
	for ( auto const & pose_id_to_indexes : incomplete_pjns_using_input_pose_ ) {
		if ( pose_id_to_indexes.second.empty() && non_deallocated_input_poses_.contains(pose_id_to_indexes.first) ) {
			final_deallocation_ids.push_back(pose_id_to_indexes.first);
		}
	}
	return final_deallocation_ids;
}

} //protocols
} //jd3
} //standard







#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::standard::PreliminaryLarvalJobTracker::save( Archive & arc ) const {
	arc( CEREAL_NVP( preliminary_job_node_inds_ ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( pjn_job_ind_begin_ ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( pjn_job_ind_end_ ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( preliminary_job_nodes_assigned_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( outstanding_job_count_ ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( pjn_max_nstruct_ ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( incomplete_pjns_using_input_pose_ ) ); // std::map<core::Size, std::set<core::Size> >
	arc( CEREAL_NVP( non_deallocated_input_poses_ ) ); // utility::vector1<core::Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::standard::PreliminaryLarvalJobTracker::load( Archive & arc ) {
	arc( preliminary_job_node_inds_ ); // utility::vector1<core::Size>
	arc( pjn_job_ind_begin_ ); // utility::vector1<core::Size>
	arc( pjn_job_ind_end_ ); // utility::vector1<core::Size>
	arc( preliminary_job_nodes_assigned_ ); // utility::vector1<_Bool>
	arc( outstanding_job_count_ ); // utility::vector1<core::Size>
	arc( pjn_max_nstruct_ ); // utility::vector1<core::Size>
	arc( incomplete_pjns_using_input_pose_ ); // std::map<core::Size, std::set<core::Size> >
	arc( non_deallocated_input_poses_ ); // utility::vector1<core::Size>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::standard::PreliminaryLarvalJobTracker );
CEREAL_REGISTER_TYPE( protocols::jd3::standard::PreliminaryLarvalJobTracker )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_standard_PreliminaryLarvalJobTracker )
#endif // SERIALIZATION
