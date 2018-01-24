// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/dag_node_managers/NodeManager.cc
/// @brief Base class for the family of JD3 Node Managers. This class is intended to offload some of the result-sorting logic from a Job Queen for a single job-dag node.
/// @detailed The base class is somewhat messy, so there are a few derived classes in protocols/jd3/DerivedNodeManagers.hh that simplify some of these interfaces by specializing for certain cases. Features include:
/// - Sorting results by some metric determined by the user (more negative values are considered "better")
/// - Partitioning results into separate bins
/// - Keeping track of jobs that have been susbmitted and the global job offset
/// - Finishing early (if result_threshold argument in the constructor is != 0) if enough results come in
/// - Determining which jobs results should be discarded and which should be kept
/// @author Jack Maguire, jack@med.unc.edu


#include <protocols/jd3/dag_node_managers/NodeManager.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.jd3.dag_node_managers.NodeManager" );

namespace protocols {
namespace jd3 {
namespace dag_node_managers {

//Constructor
NodeManager::NodeManager(
	core::Size job_offset,
	core::Size num_jobs_total,
	core::Size num_partitions,
	utility::vector1< core::Size > num_results_to_keep_for_part,
	utility::vector1< core::Size > result_threshold_per_part
) :
	job_offset_( job_offset ),
	num_jobs_total_( num_jobs_total ),

	num_results_to_keep_( 0 ),
	num_results_to_keep_for_part_( num_results_to_keep_for_part ),

	num_results_total_( 0 ),
	num_results_received_( 0 ),
	num_results_received_for_part_( num_partitions, 0 ),

	result_threshold_( ! result_threshold_per_part.empty() ),
	result_threshold_per_part_( result_threshold_per_part ),
	stopped_early_( false ),

	num_jobs_submitted_( 0 ),
	num_jobs_completed_( 0 ),

	num_partitions_( num_partitions ),
	results_to_keep_( num_partitions ),
	results_have_been_requested_at_least_once_( false )
{
	for ( core::Size part = 1; part <= num_partitions; ++part ) {
		TR.Debug << "Keeping " << num_results_to_keep_for_part_[ part ] << " results for partition " << part << std::endl;
		results_to_keep_[ part ].reserve( num_results_to_keep_for_part_[ part ] + 1 );
		num_results_to_keep_ += num_results_to_keep_for_part_[ part ];
	}

}

//Destructor
NodeManager::~NodeManager()
{}

void
NodeManager::register_result( core::Size global_job_id, core::Size local_result_id, core::Real score, core::Size partition ){
	//This does not have to be a permenant feature:
	//For now, do not touch the container of results after they have been accessed by the job queen. This is to preserve indexing.
	//Feel free to change this behavior. One option is to have two containers, 1 having all of the results that have been returned and another that has all of the results that have not been returned. When a new result comes in, sort it in the second container so that the index of the results in the first container will never change.
	if ( results_have_been_requested_at_least_once_ ) {
		job_results_that_should_be_discarded_.push_back( std::make_pair( global_job_id, local_result_id ) );
		return;
	}

	++num_results_received_;
	++num_results_received_for_part_[ partition ];

	result_elements const new_element( global_job_id, local_result_id, score );
	results_to_keep_[ partition ].insert(
		std::upper_bound( results_to_keep_[ partition ].begin(), results_to_keep_[ partition ].end(), new_element ),
		new_element
	);

	if ( results_to_keep_[ partition ].size() > num_results_to_keep_for_part_[ partition ] ) {
		result_elements const to_be_removed = results_to_keep_[ partition ].back();
		results_to_keep_[ partition ].pop_back();
		job_results_that_should_be_discarded_.push_back( std::make_pair( to_be_removed.global_job_id, to_be_removed.local_result_id ) );
	}

	if ( result_threshold_ && ready_to_finish_early() ) {
		stopped_early_ = true;
	}
}

void
NodeManager::note_job_completed( core::Size, core::Size nresults ){//first arg is global_job_id
	++num_jobs_completed_;
	//see NodeManager::register_result() comment for details on this if statement
	if ( ! results_have_been_requested_at_least_once_ ) {
		num_results_total_ += nresults;
	}
}

utility::vector1< result_elements > const &
NodeManager::results_to_keep() const {
	results_have_been_requested_at_least_once_ = true;

	all_results_to_keep_.clear();
	all_results_to_keep_.reserve( num_results_to_keep_ );
	for ( core::Size part = 1; part <= num_partitions_; ++part ) {
		all_results_to_keep_.insert( all_results_to_keep_.end(), results_to_keep_[ part ].begin(), results_to_keep_[ part ].end() );
	}
	return all_results_to_keep_;
}

jd3::JobResultID
NodeManager::get_nth_job_result_id( core::Size n ) const {
	results_have_been_requested_at_least_once_ = true;

	core::Size count = 0;
	core::Size partition = 0;
	core::Size local_index = 0;
	core::Size prev_count = 0;
	for ( core::Size ii = 1; ii <= num_partitions_; ++ii ) {
		count += num_results_to_keep_for_part_[ ii ];
		if ( n <= count ) {
			partition = ii;
			local_index = n - prev_count;
			break;
		}
		prev_count = count;
	}
	if ( local_index <= results_to_keep_[ partition ].size() ) {
		return std::make_pair( results_to_keep_[ partition ][ local_index ].global_job_id, results_to_keep_[ partition ][ local_index ].local_result_id );
	} else {
		return jd3::JobResultID(0,0);
	}
}


}
} //jd3
} //protocols
