// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/dag_node_managers/EvenlyPartitionedNodeManager.hh
/// @brief Node Manager that separates job results into different pools of equal size. Assignment of job results into pools is determined by you, the user, in the form of the "partition" argument in register_result(). If you are implementing a result threshold, this will also be evenly divided among all of the partitions
/// @detailed See here for more info: https://www.rosettacommons.org/docs/latest/development_documentation/tutorials/jd3_derived_jq/classes/node_manager
/// @author Jack Maguire, jackmaguire1444@gmail.com


#ifndef INCLUDED_protocols_jd3_dag_node_managers_EvenlyPartitionedNodeManager_HH
#define INCLUDED_protocols_jd3_dag_node_managers_EvenlyPartitionedNodeManager_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/jd3/dag_node_managers/EvenlyPartitionedNodeManager.fwd.hh>
#include <protocols/jd3/dag_node_managers/NodeManager.hh>

namespace protocols {
namespace jd3 {
namespace dag_node_managers {

//////
//UTIL
inline utility::vector1< core::Size > determine_num_for_partition(
	core::Size const num_partitions,
	core::Size const num_results_to_keep,
	bool return_empty_array_if_0
){

	if ( return_empty_array_if_0 && num_results_to_keep == 0 ) {
		return utility::vector1< core::Size > ( 0 );
	}

	utility::vector1< core::Size > num_per_partition( num_partitions, num_results_to_keep/num_partitions );
	auto const mod_result = num_results_to_keep % num_partitions;
	for ( core::Size ii = 1; ii <= mod_result; ++ii ) {
		++num_per_partition[ ii ];
	}
	return num_per_partition;
}

class EvenlyPartitionedNodeManager : public NodeManager {

public:

	//constructor
	///@param job_offset The node manager can only represent nodes where the jobids form a continuous range. That range should start at job_offset+1.
	///@param num_jobs_total The range mentioned previous should end with job_offset+num_jobs_total (inclusive)
	///@param num_results_to_keep_for_part maximum number of results you want to keep. This will divided evenly among the partitions.
	///@param num_partitions The total number of partitions you want.
	///@param result_threshold If you want to have a result threshold, define it here. This will divided evenly among the partitions.
	///We will stop submitting once every partition meets its result quota.
	EvenlyPartitionedNodeManager(
		core::Size job_offset,
		core::Size num_jobs_total,
		core::Size total_num_results_to_keep,
		core::Size num_partitions,
		core::Size result_threshold = 0,
		bool return_results_depth_first = false
	) :
		NodeManager(
		job_offset,
		num_jobs_total,
		num_partitions,
		determine_num_for_partition( num_partitions, total_num_results_to_keep, false ),
		determine_num_for_partition( num_partitions, result_threshold, true ),
		return_results_depth_first
		)
	{}

	/*
	Other useful methods:
	*/


};


}
} //jd3
} //protocols

#endif
