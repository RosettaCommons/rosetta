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

struct result_elements {
	result_elements( core::Size global_job_id_in, core::Size local_result_id_in, core::Real score_in ){
		global_job_id = global_job_id_in;
		local_result_id = local_result_id_in;
		score = score_in;
	}

	core::Size global_job_id;
	core::Size local_result_id;
	core::Real score;

	JobResultID
	job_result_id() const {
		return jd3::JobResultID( global_job_id, local_result_id );
	}

	///@brief this is for sorting
	bool operator < ( result_elements const & other ) const{
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

	//returns thing being removed
	result_elements insert( core::Size partition, result_elements const & new_guy );

	result_elements const & get_nth_element( core::Size n );
	result_elements const & get_specific_element( core::Size partition, core::Size index );

	utility::vector1< result_elements > linear_vector_of_results();

	void return_results_depth_first(){
		return_results_depth_first_ = true;
	}

	void return_results_breadth_first(){
		return_results_depth_first_ = false;
	}

	void set_return_results_depth_first( bool setting ) {
		return_results_depth_first_ = setting;
	}

private:
	utility::vector1< core::Size > n_results_to_keep_for_partition_;
	utility::vector1< core::Size > n_results_accessed_for_partition_;
	utility::vector1< utility::vector1< result_elements > > results_for_partition_;

	bool return_results_depth_first_;
};

} //dag_node_managers
} //jd3
} //protocols

#endif
