// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobGenealogist.cc
/// @brief  class method definitions for JobGenealogist
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/JobGenealogist.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/vector1.functions.hh>

namespace protocols {
namespace jd3 {

std::pair< core::Size, core::Size >
size_pair( core::Size first, core::Size second ) {
	return std::make_pair( first, second );
}

JobGenealogist::JobGenealogist() :
	num_nodes_( 0 ),
	committed_job_range_max_( 0 )
{}

JobGenealogist::JobGenealogist( JobGenealogist const & ) = default;

JobGenealogist::~JobGenealogist() = default;

JobGenealogist & JobGenealogist::operator = ( JobGenealogist const & ) = default;

void JobGenealogist::set_num_nodes( Size num_nodes )
{
	debug_assert( num_nodes_ <= num_nodes );
	if ( num_nodes_ == num_nodes ) return;

	num_nodes_ = num_nodes;
	target_njobs_for_node_.resize( num_nodes_, 0 );
	actual_njobs_for_node_.resize( num_nodes_, 0 );
	//nodes_finished_.resize( num_nodes_, false );
	target_jobid_ranges_for_node_.resize( num_nodes_, size_pair( 0, 0 ) );
	actual_jobid_ranges_for_node_.resize( num_nodes_, size_pair( 0, 0 ) );

	parent_job_groups_for_node_.resize( num_nodes_ );
	n_replicate_jobs_for_pjg_for_node_.resize( num_nodes_ );
	start_job_index_for_pjg_for_node_.resize( num_nodes_ );
	all_parental_jobs_for_node_.resize( num_nodes_ );
	pjgs_for_parental_job_for_node_.resize( num_nodes_ );

	living_jobs_w_descendants_for_node_.resize( num_nodes_ );

	last_delivered_job_for_node_.resize( num_nodes_, 0 );
	job_lookup_jobid_start_.reserve( num_nodes_ );
	job_lookup_node_index_.reserve( num_nodes_ );

	job_dag_.set_num_nodes( num_nodes_ );
}

void JobGenealogist::add_node()
{
	num_nodes_++;
	target_njobs_for_node_.resize( num_nodes_, 0 );
	actual_njobs_for_node_.resize( num_nodes_, 0 );
	//nodes_finished_.resize( num_nodes_, false );
	target_jobid_ranges_for_node_.resize( num_nodes_ );
	actual_jobid_ranges_for_node_.resize( num_nodes_ );

	parent_job_groups_for_node_.resize( num_nodes_ );
	n_replicate_jobs_for_pjg_for_node_.resize( num_nodes_ );
	start_job_index_for_pjg_for_node_.resize( num_nodes_ );
	all_parental_jobs_for_node_.resize( num_nodes_ );
	pjgs_for_parental_job_for_node_.resize( num_nodes_ );

	living_jobs_w_descendants_for_node_.resize( num_nodes_ );

	last_delivered_job_for_node_.resize( num_nodes_, 0 );
	job_lookup_jobid_start_.reserve( num_nodes_ );
	job_lookup_node_index_.reserve( num_nodes_ );

	job_dag_.add_node();
}

void JobGenealogist::set_target_num_jobs_for_node( Size node_id, Size num_jobs )
{
	debug_assert( 0 < node_id && node_id <= num_nodes_ );
	debug_assert( target_njobs_for_node_[ node_id ] == 0 );
	target_njobs_for_node_[ node_id ] = num_jobs;
	target_jobid_ranges_for_node_[ node_id ] = size_pair( committed_job_range_max_ + 1, committed_job_range_max_ + num_jobs );
	job_lookup_jobid_start_.push_back( committed_job_range_max_ + 1 );
	job_lookup_node_index_.push_back(  node_id );
	committed_job_range_max_ += num_jobs;
}

void JobGenealogist::append_parents_and_n_replicate_jobs_for_node(
	Size node_id,
	Sizes const & parent_job_index,
	Sizes const & n_replicate_jobs_for_each_parent_group
)
{
	utility::vector1< Sizes > parent_job_indices( parent_job_index.size() );
	for ( Size ii = 1; ii <= parent_job_index.size(); ++ii ) {
		parent_job_indices[ ii ].resize( 1 );
		parent_job_indices[ ii ][ 1 ] = parent_job_index[ ii ];
	}
	append_parents_and_n_replicate_jobs_for_node( node_id, parent_job_indices, n_replicate_jobs_for_each_parent_group );
}

void JobGenealogist::append_parents_and_n_replicate_jobs_for_node(
	Size node_id,
	Sizes const & parent_job_index,
	Size n_replicate_jobs_for_all_parent_groups
)
{
	Sizes all_same( parent_job_index.size(), n_replicate_jobs_for_all_parent_groups );
	append_parents_and_n_replicate_jobs_for_node( node_id, parent_job_index, all_same );
}

void JobGenealogist::append_parents_and_n_replicate_jobs_for_node(
	Size node_id,
	utility::vector1< Sizes > const & parent_job_indices,
	Sizes const & n_replicate_jobs_for_each_pjg
)
{
	Size prev_n_pjgs = parent_job_groups_for_node_[ node_id ].size();
	Size new_n_pjgs = prev_n_pjgs + parent_job_indices.size();
	parent_job_groups_for_node_[ node_id ] .resize( new_n_pjgs );
	n_replicate_jobs_for_pjg_for_node_[ node_id ].resize( new_n_pjgs );
	start_job_index_for_pjg_for_node_[ node_id ].resize( new_n_pjgs );

	Size last_parent_node = 0;
	for ( Size ii = prev_n_pjgs+1; ii <= new_n_pjgs; ++ii ) {

		for ( Size jj = 1; jj <= parent_job_indices[ ii-prev_n_pjgs ].size(); ++jj ) {
			Size const jj_job_id = parent_job_indices[ ii-prev_n_pjgs ][ jj ];
			if ( last_parent_node == 0 || ! job_is_from_node( jj_job_id, last_parent_node ) ) {
				last_parent_node = node_for_jobid( jj_job_id );
				if ( ! job_dag_.get_edge_exists( last_parent_node, node_id ) ) {
					job_dag_.add_edge( last_parent_node, node_id );
				}
			}

			if ( all_parental_jobs_for_node_[ node_id ].count( jj_job_id ) == 0 ) {
				all_parental_jobs_for_node_[ node_id ].insert( jj_job_id );
			}
			if ( pjgs_for_parental_job_for_node_[ node_id ].count( jj_job_id ) == 0 ) {
				pjgs_for_parental_job_for_node_[ node_id ][ jj_job_id ] = Sizes( 1, ii );
			} else {
				pjgs_for_parental_job_for_node_[ node_id ][ jj_job_id ].push_back( ii );
			}

			living_jobs_w_descendants_for_node_[ last_parent_node ].insert( jj_job_id );
		}

		parent_job_groups_for_node_[ node_id ][ ii ] = parent_job_indices[ ii-prev_n_pjgs ];
		n_replicate_jobs_for_pjg_for_node_[ node_id ][ ii ] = n_replicate_jobs_for_each_pjg[ ii-prev_n_pjgs ];
		actual_njobs_for_node_[ node_id ] += n_replicate_jobs_for_each_pjg[ ii-prev_n_pjgs ];
		if ( ii == 1 ) {
			start_job_index_for_pjg_for_node_[ node_id ][ ii ] = target_jobid_ranges_for_node_[ node_id ].first;
		} else {
			start_job_index_for_pjg_for_node_[ node_id ][ ii ] = start_job_index_for_pjg_for_node_[ node_id ][ ii-1 ] +
				n_replicate_jobs_for_pjg_for_node_[ node_id ][ ii-1 ];
		}
	}
	Size start_jobid = target_jobid_ranges_for_node_[ node_id ].first;
	actual_jobid_ranges_for_node_[ node_id ] = size_pair( start_jobid, start_jobid + actual_njobs_for_node_[ node_id ] - 1 );

}

void JobGenealogist::append_parents_and_n_replicate_jobs_for_node(
	Size node_id,
	utility::vector1< Sizes > const & parent_job_indices,
	Size n_replicate_jobs_for_all_parent_groups
)
{
	Sizes all_same( parent_job_indices.size(), n_replicate_jobs_for_all_parent_groups );
	append_parents_and_n_replicate_jobs_for_node( node_id, parent_job_indices, all_same );
}


void JobGenealogist::set_actual_njobs_for_node( Size node_id, Size num_jobs )
{
	debug_assert( actual_jobid_ranges_for_node_[ node_id ].first  == 0 );
	debug_assert( actual_jobid_ranges_for_node_[ node_id ].second == 0 );
	actual_njobs_for_node_[ node_id ] = num_jobs;
	Size first = target_jobid_ranges_for_node_[ node_id ].first;
	actual_jobid_ranges_for_node_[ node_id ] = size_pair( first, first+num_jobs-1 );
}

void JobGenealogist::note_job_discarded( Size job_id )
{
	//std::cout << "discarding " << job_id << std::endl;
	discarded_jobs_.insert( job_id );
	Size const node_id = node_for_jobid( job_id );

	if ( living_jobs_w_descendants_for_node_[ node_id ].count( job_id ) ) {
		living_jobs_w_descendants_for_node_[ node_id ].erase( job_id );
		for ( auto ds_node_iter = job_dag_.get_node( node_id )->outgoing_edge_list_begin(),
				ds_node_iter_end = job_dag_.get_node( node_id )->outgoing_edge_list_end();
				ds_node_iter != ds_node_iter_end; ++ds_node_iter ) {
			Size ds_node = (*ds_node_iter)->get_head_node_ind();
			all_parental_jobs_for_node_[ ds_node ].erase( job_id );
			pjgs_for_parental_job_for_node_[ ds_node ].erase( job_id );
		}
	}
}

bool JobGenealogist::job_has_been_discarded( Size job_id ) const
{
	return discarded_jobs_.member( job_id );
}

bool JobGenealogist::jobs_remain_for_node( Size node_id ) const
{
	return last_delivered_job_for_node_[ node_id ] != 0 &&
		last_delivered_job_for_node_[ node_id ] < actual_jobid_ranges_for_node_[ node_id ].second;
}

core::Size
JobGenealogist::get_next_job_for_node( Size node_id )
{
	if ( last_delivered_job_for_node_[ node_id ] == 0 ) {
		last_delivered_job_for_node_[ node_id ] = target_jobid_ranges_for_node_[ node_id ].first;
	} else if ( jobs_remain_for_node( node_id ) ) {
		++last_delivered_job_for_node_[ node_id ];
	} else {
		return 0;
	}
	return last_delivered_job_for_node_[ node_id ];
}

core::Size
JobGenealogist::get_num_nodes() const
{
	return num_nodes_;
}

core::Size JobGenealogist::get_node_target_range_begin( Size node_id ) const
{
	return target_jobid_ranges_for_node_[ node_id ].first;
}

core::Size JobGenealogist::get_node_target_range_end(   Size node_id ) const
{
	return target_jobid_ranges_for_node_[ node_id ].second;
}


core::Size JobGenealogist::get_node_actual_range_begin( Size node_id ) const
{
	return actual_jobid_ranges_for_node_[ node_id ].first;
}


core::Size JobGenealogist::get_node_actual_range_end(   Size node_id ) const
{
	return actual_jobid_ranges_for_node_[ node_id ].second;
}

bool
JobGenealogist::job_has_any_parents( Size job_id ) const
{
	Size const node_index = node_for_jobid( job_id );
	return job_from_node_has_any_parents( job_id, node_index );
}

bool
JobGenealogist::job_from_node_has_any_parents( Size job_id, Size node_id ) const
{
	if ( parent_job_groups_for_node_[ node_id ].size() == 0 ) return false;
	Size const pjg_ind = pjg_ind_for_job_from_node( job_id, node_id );
	return parent_job_groups_for_node_[ node_id ][ pjg_ind ].size() != 0;
}

core::Size
JobGenealogist::get_parent_for_job( Size job_id ) const
{
	Sizes const & parents = get_parents_for_job( job_id );

	if ( parents.size() != 1 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Tried to call get_parent_for_job (which requires that the job have only a single parent) for a job that"
			" has " + utility::to_string( parents.size() ) + " parents." );
	}

	return parents[ 1 ];
}

JobGenealogist::Sizes const &
JobGenealogist::get_parents_for_job( Size job_id ) const
{
	debug_assert( job_has_any_parents( job_id ) );
	Size node_id = node_for_jobid( job_id );
	return get_parents_for_job_for_node( job_id, node_id );
}

JobGenealogist::Sizes const &
JobGenealogist::get_parents_for_job_for_node( Size job_id, Size node_index ) const
{
	Size const pjg_ind = pjg_ind_for_job_from_node( job_id, node_index );
	return parent_job_groups_for_node_[ node_index ][ pjg_ind ];
}

JobGenealogist::Sizes
JobGenealogist::get_all_ancestors_for_job( Size job_id ) const
{
	Sizes ancestor_vector;
	Size const node_id = node_for_jobid( job_id );
	if ( !job_from_node_has_any_parents( job_id, node_id ) ) return ancestor_vector;

	std::set< core::Size > all_ancs;
	std::list< core::Size > anc_queue;

	Sizes const & parents( get_parents_for_job_for_node( job_id, node_id ));
	for ( Size ii : parents ) {
		if ( all_ancs.count( ii ) == 0 ) {
			all_ancs.insert( ii );
			anc_queue.push_back( ii );
		}
	}

	while ( !anc_queue.empty() ) {
		Size const anc_job = anc_queue.front(); anc_queue.pop_front();
		if ( ! job_has_any_parents( anc_job ) ) continue;

		Sizes const & anc_parents( get_parents_for_job_for_node( job_id, node_id ));
		for ( Size ii : anc_parents ) {
			if ( all_ancs.count( ii ) == 0 ) {
				all_ancs.insert( ii );
				anc_queue.push_back( ii );
			}
		}
	}


	ancestor_vector.resize( all_ancs.size() );
	std::copy( all_ancs.rbegin(), all_ancs.rend(), ancestor_vector.begin() );
	return ancestor_vector;
}


std::list< core::Size >
JobGenealogist::find_descendentless_jobs_backwards_from_node(
	Size finished_node_id
)
{
	// ok -- work back from the set of jobs on the finished node that have been discarded to figure out which jobs on
	// its ancestor nodes have no living descendants, and then work backwards from there, too
	std::list< Size > jobs_wo_living_descendants;

	std::list< Size > bfs_node_list;
	utility::vector1< char > node_in_bfs_list( num_nodes_, 0 );
	add_upstream_nodes_to_bfs_queue( finished_node_id, bfs_node_list, node_in_bfs_list );

	if ( bfs_node_list.empty() ) {
		return jobs_wo_living_descendants;
	}

	while ( ! bfs_node_list.empty() ) {
		Size const curr_node = bfs_node_list.front(); bfs_node_list.pop_front();
		add_upstream_nodes_to_bfs_queue( curr_node, bfs_node_list, node_in_bfs_list );

		// now iterate over all of the jobs in living_jobs_w_descendants_for_node, and
		// for each job, look at all the downstream nodes and ask if any of the jobs
		// that count this job as a parent are still alive (undiscarded)
		std::list< Size > new_jobs_wo_living_descendants;
		for ( auto job_id : living_jobs_w_descendants_for_node_[ curr_node ] ) {
			bool found_living_descendant = false;
			//bool keep_going = true;
			for ( auto ds_node_iter = job_dag_.get_node( curr_node )->outgoing_edge_list_begin(),
					ds_node_iter_end = job_dag_.get_node( curr_node )->outgoing_edge_list_end();
					ds_node_iter != ds_node_iter_end; ++ds_node_iter ) {
				Size ds_node = (*ds_node_iter)->get_head_node_ind();
				if ( ! all_parental_jobs_for_node_[ ds_node ].count( job_id ) ) continue;
				for ( Size ii : pjgs_for_parental_job_for_node_[ ds_node ][ job_id ] ) {
					for ( Size jj = start_job_index_for_pjg_for_node_[ ds_node ][ ii ];
							jj <= start_job_index_for_pjg_for_node_[ ds_node ][ ii ] +
							n_replicate_jobs_for_pjg_for_node_[ ds_node ][ ii ] - 1; ++jj ) {
						if ( ! discarded_jobs_.member( jj ) ) {
							found_living_descendant = true;
							break;
						}
					}
					if ( found_living_descendant ) break;
				}
				if ( found_living_descendant ) break;
			}
			if ( ! found_living_descendant ) { new_jobs_wo_living_descendants.push_back( job_id ); }
		}

		// Now iterate across all the jobs we just collected that have no living descendants
		// and erase record of them in the downstream nodes so we don't again spend time
		// trying to decide if they have any living descendants.
		if ( ! new_jobs_wo_living_descendants.empty() ) {
			for ( auto ds_node_iter = job_dag_.get_node( curr_node )->outgoing_edge_list_begin(),
					ds_node_iter_end = job_dag_.get_node( curr_node )->outgoing_edge_list_end();
					ds_node_iter != ds_node_iter_end; ++ds_node_iter ) {
				Size ds_node = (*ds_node_iter)->get_head_node_ind();
				for ( Size ii : new_jobs_wo_living_descendants ) {
					all_parental_jobs_for_node_[ ds_node ].erase( ii );
					pjgs_for_parental_job_for_node_[ ds_node ].erase( ii );
				}
			}
		}

		// now append the newly discovered jobs to the growing list
		jobs_wo_living_descendants.splice( jobs_wo_living_descendants.end(), new_jobs_wo_living_descendants );
	}

	return jobs_wo_living_descendants;
}

core::Size
JobGenealogist::node_for_jobid( Size job_id ) const
{
	Size const ind = utility::binary_search_ranges( job_lookup_jobid_start_, job_id );
	Size const node_index = job_lookup_node_index_[ ind ];
	debug_assert( target_jobid_ranges_for_node_[ node_index ].first  <= job_id );
	debug_assert( target_jobid_ranges_for_node_[ node_index ].second >= job_id );
	return node_index;
}

core::Size
JobGenealogist::replicate_id_for_jobid( Size job_id ) const
{
	std::pair< core::Size, core::Size > node_n_repid = node_and_replicate_id_for_jobid( job_id );
	return node_n_repid.second;
}

std::pair< core::Size, core::Size >
JobGenealogist::node_and_replicate_id_for_jobid( Size job_id ) const
{
	Size const node_index = node_for_jobid( job_id );
	Size const replicate_id = replicate_id_for_jobid_from_node( job_id, node_index );
	return std::make_pair( node_index, replicate_id );
}

core::Size
JobGenealogist::replicate_id_for_jobid_from_node( Size job_id, Size node_index ) const
{
	debug_assert( job_is_from_node( job_id, node_index ) );
	if ( parent_job_groups_for_node_[ node_index ].size() != 0 ) {
		Size const pjg_ind = pjg_ind_for_job_from_node( job_id, node_index );
		return job_id - start_job_index_for_pjg_for_node_[ node_index ][ pjg_ind ] + 1;
	} else {
		return job_id - target_jobid_ranges_for_node_[ node_index ].first + 1;
	}
}

core::Size
JobGenealogist::pjg_ind_for_job_from_node( Size job_id, Size node_id ) const
{
	debug_assert( parent_job_groups_for_node_[ node_id ].size() != 0 );
	Sizes const & start_job_indices( start_job_index_for_pjg_for_node_[ node_id ] );
	return utility::binary_search_ranges( start_job_indices, job_id );
}


bool JobGenealogist::job_is_from_node( Size job_id, Size node_id ) const
{
	return actual_jobid_ranges_for_node_[ node_id ].first <= job_id && job_id < actual_jobid_ranges_for_node_[ node_id ].second;
}

void
JobGenealogist::add_upstream_nodes_to_bfs_queue(
	Size const target_node,
	std::list< Size > & bfs_node_list,
	utility::vector1< char > & node_in_bfs_list
) const
{
	for ( auto upstream_iter = job_dag_.get_node( target_node )->const_incoming_edge_list_begin(),
			upstream_iter_end = job_dag_.get_node( target_node )->const_incoming_edge_list_end();
			upstream_iter != upstream_iter_end; ++upstream_iter ) {
		Size upstream_index = (*upstream_iter)->get_tail_node_ind();
		if ( ! node_in_bfs_list[ upstream_index ] ) {
			bfs_node_list.push_back( upstream_index );
			node_in_bfs_list[ upstream_index ] = 1;
		}
	}

}

} // namespace jd3
} // namespace protocols

