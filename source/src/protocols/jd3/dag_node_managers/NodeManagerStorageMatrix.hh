// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/dag_node_managers/NodeManagerStorageMatrix.hh
/// @author Jack Maguire, jackmaguire1444@gmail.com

/*
One of the options for this class involves the order in which results are returned
accross partitions. To get a better understanding, consider the following example:

..............Results

Partition 1: A B C D

Partition 2: E F G H

Partition 3: I J K L

Breadth-First Order: A E I B F J C G K D H L

Depth-First Order:   A B C D E F G H I J K L

*/

#ifndef INCLUDED_protocols_jd3_dag_node_managers_NodeManagerStorageMatrix_HH
#define INCLUDED_protocols_jd3_dag_node_managers_NodeManagerStorageMatrix_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/jd3/dag_node_managers/NodeManagerStorageMatrix.fwd.hh>
#include <protocols/jd3/CompletedJobOutput.fwd.hh>
//#include <protocols/jd3/dag_node_managers/NodeManager.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace jd3 {
namespace dag_node_managers {

struct ResultElements {
	ResultElements(
		core::Size global_job_id_in,
		core::Size local_result_id_in,
		core::Real score_in,
		uint64_t token_in = 0
	) :
		global_job_id( global_job_id_in ),
		local_result_id( local_result_id_in ),
		score( score_in ),
		token( token_in )
	{}

	core::Size global_job_id;
	core::Size local_result_id;
	core::Real score;
	uint64_t token;

	JobResultID
	job_result_id() const {
		return jd3::JobResultID( global_job_id, local_result_id );
	}

	///@brief this is for sorting
	bool operator < ( ResultElements const & other ) const{
		return score < other.score;
	}
};


class NodeManagerStorageMatrix {

public:
	NodeManagerStorageMatrix(
		utility::vector1< core::Size > n_results_to_keep_for_partition,//by-value because we std::move the copy
		bool return_results_depth_first = false
	);

	virtual ~NodeManagerStorageMatrix();

	//returns the element being displaced/removed if applicable
	ResultElements insert( core::Size partition, ResultElements const & new_guy );

	ResultElements const & get_nth_element( core::Size n );
	ResultElements const & get_specific_element( core::Size partition, core::Size index_within_partition );

	///@brief Freezes the dynamic arrays and appends them onto each other (using the depth-first pattern)
	utility::vector1< ResultElements >
	linear_vector_of_results();

	///@brief Consider example at the top of this file to learn more about this option
	void set_return_results_depth_first( bool setting ) {
		return_results_depth_first_ = setting;
	}

	void return_results_depth_first(){
		return_results_depth_first_ = true;
	}

	void return_results_breadth_first(){
		return_results_depth_first_ = false;
	}

	///@brief Perhaps you have a diversity requirement and do not want too many results with the same token
	/// (token can represent anything you want - as long as it can be stored as a uint64_t), this setting is for you.
	void set_max_num_results_with_same_token_per_partition( core::Size setting ) {
		max_num_results_with_same_token_per_partition_ = setting;
	}

	core::Size num_partitions() const {
		core::Size const n_part = n_results_to_keep_for_partition_.size();
		debug_assert( n_results_accessed_for_partition_.size() == n_part );
		debug_assert( results_for_partition_.size() == n_part );
		return n_part;
	}

	void clear(){
		//n_results_to_keep_for_partition_.clear();
		//n_results_accessed_for_partition_.clear();
		results_for_partition_.clear();
	}

private:
	utility::vector1< core::Size > n_results_to_keep_for_partition_;
	utility::vector1< core::Size > n_results_accessed_for_partition_;
	utility::vector1< utility::vector1< ResultElements > > results_for_partition_;

	bool return_results_depth_first_;

	core::Size max_num_results_with_same_token_per_partition_;
};

} //dag_node_managers
} //jd3
} //protocols

#endif
