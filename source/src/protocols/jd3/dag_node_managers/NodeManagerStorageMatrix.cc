// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/dag_node_managers/NodeManagerStorageMatrix.cc
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/jd3/dag_node_managers/NodeManagerStorageMatrix.hh>
#include <algorithm>

namespace protocols {
namespace jd3 {
namespace dag_node_managers {

NodeManagerStorageMatrix::NodeManagerStorageMatrix(
	utility::vector1< core::Size > n_results_to_keep_for_partition,//by-value because we std::move the copy
	bool return_results_depth_first
) :
	n_results_to_keep_for_partition_( std::move( n_results_to_keep_for_partition ) ),
	n_results_accessed_for_partition_( n_results_to_keep_for_partition_.size(), 0 ),
	results_for_partition_( n_results_to_keep_for_partition_.size() ),
	return_results_depth_first_( return_results_depth_first ),
	max_num_results_with_same_token_per_partition_( 0 )
{
	core::Size const n_partitions = n_results_to_keep_for_partition_.size();
	for ( core::Size part = 1; part <= n_partitions; ++part ) {
		//Reserve an extra element for temporary overflow during insertion but before removal
		results_for_partition_[ part ].reserve( n_results_to_keep_for_partition_[ part ] + 1 );
	}
}

NodeManagerStorageMatrix::~NodeManagerStorageMatrix(){}

core::Size
get_index_of_last_element_with_token (
	utility::vector1< ResultElements > const & vec,
	core::Size token
){
	for ( auto ii = vec.size(); ii > 0; --ii ) {
		if ( vec[ ii ].token == token ) return ii;
	}
	utility_exit_with_message( "token: " + std::to_string(token) + " is not in vector." );
	return 0;
}

ResultElements
NodeManagerStorageMatrix::insert( core::Size const partition, ResultElements const & new_guy ){

	utility::vector1< ResultElements > & partition_vec = results_for_partition_[ partition ];
	ResultElements doomed_guy( 0, 0, 0 );

	//First determine if this exceeds our token limit
	if ( max_num_results_with_same_token_per_partition_ != 0 ) {
		auto pred = [ & new_guy ]( ResultElements const & i ) -> bool { return i.token == new_guy.token; };
		core::Size const num_elements_with_this_token =
			std::count_if( partition_vec.begin(), partition_vec.end(), pred );

		if ( num_elements_with_this_token >= max_num_results_with_same_token_per_partition_ ) {
			//compare new_guy to the lowest rank result with the same element
			unsigned int const index_of_last_token = get_index_of_last_element_with_token( partition_vec, new_guy.token );
			if ( index_of_last_token <= n_results_accessed_for_partition_[ partition ] ) {
				//region is frozen, can not add new guy
				return new_guy;
			} if ( new_guy.score >= partition_vec[ index_of_last_token ].score ) {
				//new guy is not good enough
				return new_guy;
			} else {
				//new guy is good enough. Delete
				doomed_guy = partition_vec[ index_of_last_token ];
				partition_vec.erase( partition_vec.begin() + index_of_last_token - 1 );
			}
		}
	}

	utility::vector1< ResultElements >::iterator first_available_position =
		partition_vec.begin() + n_results_accessed_for_partition_[ partition ];

	partition_vec.insert(
		std::upper_bound( first_available_position, partition_vec.end(), new_guy ),
		new_guy
	);

	if ( partition_vec.size() > n_results_to_keep_for_partition_[ partition ] ) {
		ResultElements const return_val = partition_vec.back();
		partition_vec.pop_back();
		debug_assert( ! ( doomed_guy.global_job_id || doomed_guy.local_result_id ) );
		return return_val;
	} else {
		return doomed_guy;
	}

}

ResultElements const &
NodeManagerStorageMatrix::get_nth_element( core::Size n ){

	core::Size const n_partitions = n_results_accessed_for_partition_.size();

	if ( return_results_depth_first_ ) {
		for ( core::Size part = 1; part <= n_partitions; ++part ) {
			if ( n_results_to_keep_for_partition_[ part ] < n ) {
				n -= n_results_to_keep_for_partition_[ part ];
			} else {
				if ( n > n_results_accessed_for_partition_[ part ] ) {
					//update freeze pointer
					n_results_accessed_for_partition_[ part ] = n;
				}

				while ( n > results_for_partition_[ part ].size() ) {
					//pad with blanks so that more results can come in the future
					results_for_partition_[ part ].emplace_back( 0, 0, 0 );
				}

				return results_for_partition_[ part ][ n ];
			}
		}
		utility_exit_with_message( "dead_code" );
		//return ResultElements( 0, 0, 0 );
	} else {
		core::Size const minimum_num_results_to_keep =
			* std::min_element( n_results_to_keep_for_partition_.begin(), n_results_to_keep_for_partition_.end() );
		core::Size const rectangle_area = minimum_num_results_to_keep * n_partitions;

		if ( n <= rectangle_area ) {
			auto const n_minus_1 = n - 1;
			auto const zero_indexed_partition = n_minus_1 % n_partitions;
			auto const zero_indexed_element = n_minus_1 / n_partitions;

			auto const part = zero_indexed_partition + 1;
			auto const element = zero_indexed_element + 1;

			if ( element > n_results_accessed_for_partition_[ part ] ) {
				//update freeze pointer
				n_results_accessed_for_partition_[ part ] = element;
			}

			while ( element > results_for_partition_[ part ].size() ) {
				//pad with blanks so that more results can come in the future
				results_for_partition_[ part ].emplace_back( 0, 0, 0 );
			}

			return results_for_partition_[ part ][ element ];

		} else {
			n -= rectangle_area;

			core::Size part = 0;
			core::Size element = minimum_num_results_to_keep + 1;

			while ( n > 0 ) {
				++part;
				if ( part > n_partitions ) {
					part = 1;
					++element;
				}
				if ( element <= n_results_to_keep_for_partition_[ part ] ) {
					--n;
				}
			}

			runtime_assert( part );

			if ( element > n_results_accessed_for_partition_[ part ] ) {
				//update freeze pointer
				n_results_accessed_for_partition_[ part ] = element;
			}

			while ( element > results_for_partition_[ part ].size() ) {
				//pad with blanks so that more results can come in the future
				results_for_partition_[ part ].emplace_back( 0, 0, 0 );
			}

			return results_for_partition_[ part ][ element ];
		}
	}

}

ResultElements const &
NodeManagerStorageMatrix::get_specific_element( core::Size partition, core::Size index ){
	while ( index > results_for_partition_[ partition ].size() ) {
		//pad with blanks so that more results can come in the future
		results_for_partition_[ partition ].emplace_back( 0, 0, 0 );
	}


	if ( index > n_results_accessed_for_partition_[ partition ] ) {
		n_results_accessed_for_partition_[ partition ] = index;
	}

	return results_for_partition_[ partition ][ index ];
}

utility::vector1< ResultElements >
NodeManagerStorageMatrix::linear_vector_of_results() {
	core::Size nresults = 0;
	for ( utility::vector1< ResultElements > const & subvec : results_for_partition_ ) {
		nresults += subvec.size();
	}
	utility::vector1< ResultElements > all_results;
	all_results.reserve( nresults );

	core::Size const n_partitions = n_results_to_keep_for_partition_.size();
	for ( core::Size part = 1; part <= n_partitions; ++part ) {
		utility::vector1< ResultElements > const & subvec = results_for_partition_[ part ];
		n_results_accessed_for_partition_[ part ] = subvec.size();
		all_results.insert( all_results.end(), subvec.begin(), subvec.end() );
	}

	runtime_assert( all_results.size() == nresults );
	return all_results;
}

} //dag_node_managers
} //jd3
} //protocols
