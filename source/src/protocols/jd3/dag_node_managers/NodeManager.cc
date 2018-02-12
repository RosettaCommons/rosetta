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
/// @author Jack Maguire, jackmaguire1444@gmail.com


#include <protocols/jd3/dag_node_managers/NodeManager.hh>
#include <protocols/jd3/dag_node_managers/NodeManagerStorageMatrix.hh>

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
	utility::vector1< core::Size > result_threshold_per_part,
	bool return_results_depth_first
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
	results_to_keep_( num_results_to_keep_for_part_, return_results_depth_first )
{
	for ( core::Size part = 1; part <= num_partitions; ++part ) {
		TR.Debug << "Keeping " << num_results_to_keep_for_part_[ part ] << " results for partition " << part << std::endl;
		num_results_to_keep_ += num_results_to_keep_for_part_[ part ];
	}

}

//Destructor
NodeManager::~NodeManager()
{}

void
NodeManager::register_result( core::Size global_job_id, core::Size local_result_id, core::Real score, core::Size partition ){

	++num_results_received_;
	++num_results_received_for_part_[ partition ];

	result_elements const new_element( global_job_id, local_result_id, score );
	result_elements const elem_to_delete = results_to_keep_.insert( partition, new_element );

	if ( elem_to_delete.global_job_id || elem_to_delete.local_result_id ) {
		//Do not consider {0,0}
		job_results_that_should_be_discarded_.emplace_back( elem_to_delete.global_job_id, elem_to_delete.local_result_id );
	}

	if ( result_threshold_ && ready_to_finish_early() ) {
		stopped_early_ = true;
	}
}

}
} //jd3
} //protocols
