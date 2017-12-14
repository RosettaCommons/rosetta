// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pack_daemon/MultistateFitnessFunction.cc
/// @brief  Implementation of the class MultistateFitnessFunction
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/pack_daemon/MultistateFitnessFunction.hh>

// Package headers
#include <protocols/pack_daemon/MultistateAggregateFunction.hh>
#include <protocols/pack_daemon/PackDaemon.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/assert.hh>
//#include <utility/heap.hh>
#include <utility/mpi_util.hh>
#include <utility/string_util.hh>

// C++ headers
#include <iostream>

#include <core/import_pose/import_pose.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


static basic::Tracer TR( "protocols.pack_daemon.MultistateFitnessFunction" );

namespace protocols {
namespace pack_daemon {


/// A struct for compairing entities and state energies for use in the
/// STL heap functions
struct EntityHistoryLT
{
public:
	using EntityAndScore = TopEntitySet::EntityAndScore;

public:
	bool operator () ( EntityAndScore const & a, EntityAndScore const & b ) {
		return a.first->fitness() < b.first->fitness();
	}
};

TopEntitySet::TopEntitySet() :
	desired_entity_history_size_( 0 ),
	n_tied_for_worst_( 0 )
{
}

core::Size
TopEntitySet::size() const
{
	return top_entities_.size();
}

TopEntitySet::EntityAndScore const &
TopEntitySet::operator[] ( core::Size index ) const {
	return top_entities_[ index ];
}

TopEntitySet::EntityAndScore &
TopEntitySet::operator[] ( core::Size index ) {
	return top_entities_[ index ];
}

void
TopEntitySet::desired_entity_history_size( core::Size setting )
{
	desired_entity_history_size_ = setting;
	top_entities_.clear();
	top_entities_.reserve( 2 * desired_entity_history_size_ );
	n_tied_for_worst_ = 0;
}

void
TopEntitySet::clear()
{
	top_entities_.clear();
	n_tied_for_worst_ = 0;
}

core::Size
TopEntitySet::desired_entity_history_size() const
{
	return desired_entity_history_size_;
}

std::list< TopEntitySet::EntityOP >
TopEntitySet::update_entity_history(
	Entity const & ent,
	StateEnergiesAndNPDs const & seanpds,
	bool & added_new_entity
)
{
	added_new_entity = false;
	std::list< EntityOP > removed_entities;

	if ( top_entities_.size() < desired_entity_history_size_ ) {
		// 1. handle the cases when the heap is not yet full

		if ( top_entities_.size() == 0 ) {
			// a. first element added to the heap
			n_tied_for_worst_ = 1;
		} else if ( ent.fitness() == top_entities_.front().first->fitness() ) {
			// b. tied with the worst element in the heap
			++n_tied_for_worst_;
		} else if ( ent.fitness() > top_entities_.front().first->fitness() ) {
			// c. worse than the worst element in the heap -- a new low!
			n_tied_for_worst_ = 1;
		}

		top_entities_.push_back( std::make_pair( EntityOP( new Entity( ent ) ), seanpds ) );
		std::push_heap( top_entities_.begin(), top_entities_.end(), EntityHistoryLT() );
		added_new_entity = true;

	} else if ( ent.fitness() <= top_entities_.front().first->fitness() ) {
		// 2. handle the cases when the heap is full
		if ( top_entities_.front().first->fitness() == ent.fitness() ) {
			// a. Tie between the new entity and the worst entity in the heap
			++n_tied_for_worst_;
		} else {
			// b. Not a tie -- maybe we can remove some of the entities that are tied for worst
			debug_assert( n_tied_for_worst_ <= top_entities_.size() );
			if ( top_entities_.size() + 1 - n_tied_for_worst_ >= desired_entity_history_size_ ) {
				// start popping entities from the heap
				ASSERT_ONLY( core::Real const old_worst_fitness = top_entities_.front().first->fitness(););
				for ( Size ii = 1; ii <= n_tied_for_worst_; ++ii ) {
					debug_assert( top_entities_.front().first->fitness() == old_worst_fitness );
					removed_entities.push_back( top_entities_.front().first );
					std::pop_heap( top_entities_.begin(), top_entities_.end(), EntityHistoryLT() );
					top_entities_.pop_back();
				}

				/// now, if the new entity has the same fitness as the current worst entity, determine how many entities are tied for worst.
				if ( top_entities_.size() != 0 && ent.fitness() == top_entities_.front().first->fitness() ) {
					utility::vector1< EntityAndScore > worst_entities; worst_entities.reserve( top_entities_.size() );
					core::Real const new_worst_fitness = top_entities_.front().first->fitness();
					n_tied_for_worst_ = 1; // count from 1, because the entity we're about to add will also be among the set of the worst
					while ( top_entities_.size() != 0 && top_entities_.front().first->fitness() == new_worst_fitness ) {
						worst_entities.push_back( top_entities_.front() );
						std::pop_heap( top_entities_.begin(), top_entities_.end(), EntityHistoryLT() );
						top_entities_.pop_back();
						++n_tied_for_worst_;
					}
					/// ok -- we know how many entities are tied for worst.  put these entities back
					/// into the heap
					for ( Size ii = 1; ii <= worst_entities.size(); ++ii ) {
						top_entities_.push_back( worst_entities[ ii ] );
						std::push_heap( top_entities_.begin(), top_entities_.end(), EntityHistoryLT() );
					}
				} else {
					/// we just emptied out the history; the element we're about to add
					/// will be the worst element in the history -- so save the fact that the number
					/// of entities tied for last place is now 1.
					n_tied_for_worst_ = 1;
				}
			}

		}
		top_entities_.push_back( std::make_pair( EntityOP( new Entity( ent ) ), seanpds ) );
		std::push_heap( top_entities_.begin(), top_entities_.end(), EntityHistoryLT() );
		added_new_entity = true;
	}

	return removed_entities;
}

core::Size
TopEntitySet::index_of_entity( Entity const & ent ) const
{
	for ( core::Size ii = 1; ii <= top_entities_.size(); ++ii ) {
		if ( *top_entities_[ ii ].first == ent ) {
			return ii;
		}
	}
	return 0;
}


TopEntitySet::EntityAndScore
TopEntitySet::pop()
{
	if ( top_entities_.size() == 0 ) {
		EntityAndScore empty;
		return empty;
	}
	EntityAndScore worst = top_entities_.front();
	std::pop_heap( top_entities_.begin(), top_entities_.end(), EntityHistoryLT() );
	return worst;
}


MultistateFitnessFunction::MultistateFitnessFunction()
{
	set_history_size( 1 ); // at the very least, store the optimal solution!
}
MultistateFitnessFunction::~MultistateFitnessFunction() = default;

core::Real MultistateFitnessFunction::evaluate( Entity & entity )
{
	clock_t start_time = clock();

	std::fill( state_energies_.begin(), state_energies_.end(), 0.0 );
	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << "Evaluating entity " << entity << std::endl;
	}

	compute_state_energies( entity );

	core::Real score = compute_aggregate_score( entity );
	entity.set_fitness( score );
	update_entity_history( entity );

	clock_t stop_time = clock();
	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << "evaluate() took " << ((double) stop_time - start_time ) / CLOCKS_PER_SEC << " seconds" << std::endl;
	}
	return score;
}

void MultistateFitnessFunction::daemon_set( DaemonSetOP ds )
{
	daemon_set_ = ds;
	state_energies_.resize( daemon_set_->ndaemons() );
	npd_properties_.resize( daemon_set_->n_npd_properties() );
	std::fill( state_energies_.begin(), state_energies_.end(), 0.0 );
	std::fill( npd_properties_.begin(), npd_properties_.end(), 0.0 );
}

void MultistateFitnessFunction::aggregate_function( MultistateAggregateFunctionOP func )
{
	aggregate_ = func;
}

DaemonSetCOP MultistateFitnessFunction::daemon_set() const
{
	return daemon_set_;
}

utility::vector1< core::Real > const & MultistateFitnessFunction::state_energies() const
{
	return state_energies_;
}

utility::vector1< core::Real > const & MultistateFitnessFunction::npd_properties() const
{
	return npd_properties_;
}


MultistateAggregateFunctionCOP MultistateFitnessFunction::aggregate_function() const
{
	return aggregate_;
}

void MultistateFitnessFunction::set_history_size( core::Size history_size )
{
	top_entities_.desired_entity_history_size( history_size );
}

void MultistateFitnessFunction::clear_history()
{
	top_entities_.clear();
}

std::list< std::pair< MultistateFitnessFunction::Size, MultistateFitnessFunction::PoseOP > >
MultistateFitnessFunction::recover_relevant_poses_for_entity( Entity const & ent )
{
	Size entity_index = top_entities_.index_of_entity( ent );
	if ( entity_index == 0 ) {
		std::cerr << "Failed to find desired top entity: " << ent << std::endl;
		std::cerr << "Top entities:" << std::endl;
		for ( Size ii = 1; ii <= top_entities_.size(); ++ii ) {
			std::cerr << ii << " " << *top_entities_[ ii ].first << std::endl;
		}
		utility_exit_with_message( "Failed to find desired top entity" );
	}
	MultistateAggregateFunction::StateIndices relevant_states =
		aggregate_->select_relevant_states(
		top_entities_[ entity_index ].second.first,
		top_entities_[ entity_index ].second.second,
		ent );
	return recover_poses_from_states( *top_entities_[ entity_index ].first, relevant_states );
}


void MultistateFitnessFunction::compute_state_energies( Entity const & entity )
{
	//std::fill( state_energies_.begin(), state_energies_.end(), 1234 );
	//std::fill( npd_properties_.begin(), npd_properties_.end(), 1234 );
	DaemonSet::StateEsAndNPDs energies =
		daemon_set_->compute_energy_for_assignment( entity );
	for ( DaemonSet::SizeRealPairs::const_iterator
			iter = energies.first.begin(), iter_end = energies.first.end();
			iter != iter_end; ++iter ) {
		state_energies_[ iter->first ] = iter->second;
	}
	for ( DaemonSet::SizeRealPairs::const_iterator
			iter = energies.second.begin(), iter_end = energies.second.end();
			iter != iter_end; ++iter ) {
		npd_properties_[ iter->first ] = iter->second;
	}

	//for ( Size ii = 1; ii <= state_energies_.size(); ++ii ) {
	// if ( state_energies_[ ii ] == 1234 ) {
	//  std::cout << "Weird: state_energies_[ " << ii << " ] was not eevaluated" << std::endl;
	// }
	//}
	//for ( Size ii = 1; ii <= npd_properties_.size(); ++ii ) {
	// if ( npd_properties_[ ii ] == 1234 ) {
	//  std::cout << "Weird: npd_properties_[ " << ii << " ] was not eevaluated" << std::endl;
	// }
	//}
}

core::Real MultistateFitnessFunction::compute_aggregate_score( Entity const & entity )
{
	return aggregate_->evaluate( state_energies_, npd_properties_, entity );
}

void MultistateFitnessFunction::instruct_daemons_to_keep_last_entity()
{
	daemon_set_->mark_last_entity_as_important();
}

void MultistateFitnessFunction::instruct_daemons_to_drop_entity( Entity const & entity )
{
	daemon_set_->mark_entity_as_unimportant( entity );
}

std::list< std::pair< MultistateFitnessFunction::Size, MultistateFitnessFunction::PoseOP > >
MultistateFitnessFunction::recover_poses_from_states(
	Entity const & ent,
	utility::vector1< core::Size > const & which_states
)
{
	return daemon_set_->retrieve_relevant_poses_for_entity( ent, which_states );
}


utility::vector1< core::Real > &
MultistateFitnessFunction::state_energies() {
	return state_energies_;
}

utility::vector1< core::Real > &
MultistateFitnessFunction::npd_properties()
{
	return npd_properties_;
}

DaemonSetOP MultistateFitnessFunction::daemon_set()
{
	return daemon_set_;
}

MultistateAggregateFunctionOP MultistateFitnessFunction::aggregate_function()
{
	return aggregate_;
}


void
MultistateFitnessFunction::update_entity_history( Entity const & ent )
{
	StateEnergiesAndNPDs seanpd = make_pair( state_energies_, npd_properties_ );
	bool added_entity( false );
	std::list< genetic_algorithm::EntityOP > discarded = top_entities_.update_entity_history( ent, seanpd, added_entity );

	if ( added_entity ) instruct_daemons_to_keep_last_entity();

	for ( auto & iter : discarded ) {
		instruct_daemons_to_drop_entity( *iter );
	}

}

MPIMultistateFitnessFunction::MPIMultistateFitnessFunction() :
	MPI_nprocs_( 0 )
{
#ifdef USEMPI
	MPI_nprocs_ = static_cast< Size > ( utility::mpi_nprocs() );
#endif
#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	utilization_by_node_.resize( MPI_nprocs_ );
	packing_percentage_.resize( MPI_nprocs_ );
	npd_percentage_.resize( MPI_nprocs_ );
#endif
}

MPIMultistateFitnessFunction::~MPIMultistateFitnessFunction() = default;

void MPIMultistateFitnessFunction::set_num_pack_daemons( Size n_daemons )
{
	state_energies().resize( n_daemons );
}

void MPIMultistateFitnessFunction::set_num_npd_properties( Size n_npd_properties )
{
	npd_properties().resize( n_npd_properties );
}

void MPIMultistateFitnessFunction::send_spin_down_signal()
{
	for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
		utility::send_integer_to_node( ii, spin_down );
	}
}

#ifdef APL_MEASURE_MSD_LOAD_BALANCE
void MPIMultistateFitnessFunction::print_load_balance_statistics( std::ostream & ostr ) const
{
	for ( Size ii = 0; ii < MPI_nprocs_; ++ii ) {
		Real ut_sum = 0; Real packing_sum = 0; Real npd_sum = 0; Size count = 0;
		for ( std::list< core::Real >::const_iterator
				uiter = utilization_by_node_[ ii ].begin(),
				uiter_end = utilization_by_node_[ ii ].end(),
				piter = packing_percentage_[ ii ].begin(),
				//piter_end = packing_percentage_[ ii ].end(),
				npditer = npd_percentage_[ ii ].begin() //,
				//npditer_end = npd_percentage_[ ii ].end
				; uiter != uiter_end; ++uiter, ++piter, ++npditer ) {
			++count;
			ut_sum += *uiter;
			packing_sum += *piter;
			npd_sum += *npditer;
		}
		ut_sum /= count;
		packing_sum /= count;
		npd_sum /= count;
		ostr << "Node " << ii << " average utilization: " << ut_sum << " packing: " << packing_sum << " npd: " << npd_sum << std::endl;
	}
}

void MPIMultistateFitnessFunction::reset_load_balance_statistics()
{
	for ( Size ii = 0; ii < MPI_nprocs_; ++ii ) {
		utilization_by_node_[ ii ].clear();
	}
}
#endif


void MPIMultistateFitnessFunction::compute_state_energies( Entity const & entity )
{
	using namespace utility;
#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	std::clock_t start_time = clock();
	utility::vector0< core::Real > times( MPI_nprocs_ );
	utility::vector0< core::Real > packing_times( MPI_nprocs_ );
	utility::vector0< core::Real > npd_times( MPI_nprocs_ );
#endif

	for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
		send_integer_to_node( ii, evaluate_entity );
	}
	broadcast_entity_string( entity );
	std::fill( state_energies().begin(), state_energies().end(), core::Real( 0.0 ) );
	if ( daemon_set() ) {

		DaemonSet::StateEsAndNPDs energies = daemon_set()->
			compute_energy_for_assignment( entity );
		for ( std::list< std::pair< core::Size, core::Real > >::const_iterator
				iter = energies.first.begin(), iter_end = energies.first.end();
				iter != iter_end; ++iter ) {
			state_energies()[ iter->first ] = iter->second;
		}
		for ( std::list< std::pair< core::Size, core::Real > >::const_iterator
				iter = energies.second.begin(), iter_end = energies.second.end();
				iter != iter_end; ++iter ) {
			npd_properties()[ iter->first ] = iter->second;
		}
	}
#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	std::clock_t node0_stop_time = clock();
	times[ 0 ] = ((double) node0_stop_time - start_time ) / CLOCKS_PER_SEC;
	packing_times[ 0 ] = daemon_set()->last_packing_runtime();
	npd_times[ 0 ] = daemon_set()->last_npd_runtime();
#endif

	for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
		//TR << "receiving n_results from node " << ii << std::endl;
		int n_results = receive_integer_from_node( ii );
		//TR << "received " << n_results << " from node " << ii << std::endl;
		for ( int jj = 1; jj <= n_results; ++jj ) {
			int which_pack_daemon  = receive_integer_from_node( ii );
			Real energy_for_daemon = receive_double_from_node(  ii );
			//std::cout << "Received energy of " << energy_for_daemon << " for daemon #" << which_pack_daemon << " on node " << ii << " from a total of " << state_energies().size() << " PackDaemons" << std::endl;
			state_energies()[ which_pack_daemon ] = energy_for_daemon;
		}
		//TR << "receiving n npd properties from node " << ii << std::endl;
		int n_npd_properties = receive_integer_from_node( ii );
		//TR << "received " << n_npd_properties << " n_npd_properties from node " << ii << std::endl;
		for ( int jj = 1; jj <= n_npd_properties; ++jj ) {
			int which_npd_prop = receive_integer_from_node( ii );
			Real npd_property  = receive_double_from_node(  ii );
			npd_properties()[ which_npd_prop ] = npd_property;
			//TR << "property " << which_npd_prop << ": " << npd_property << " from node " << ii << std::endl;
		}
#ifdef APL_MEASURE_MSD_LOAD_BALANCE
		times[ ii ] = utility::receive_double_from_node( ii );
		packing_times[ ii ] = utility::receive_double_from_node( ii );
		npd_times[ ii ] = utility::receive_double_from_node( ii );
#endif
	}


#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	std::clock_t final_stop_time = clock();
	Real runtime = ((double) final_stop_time - start_time ) / CLOCKS_PER_SEC;
	if (runtime != 0.0 ) {
		for ( Size ii = 0; ii < MPI_nprocs_; ++ii ) {
			utilization_by_node_[ ii ].push_back( times[ ii ] / runtime );
			packing_percentage_[ ii ].push_back( packing_times[ ii ] / runtime );
			npd_percentage_[ ii ].push_back( npd_times[ ii ] / runtime );
		}
	} else {
		for ( Size ii = 0; ii < MPI_nprocs_; ++ii ) {
			utilization_by_node_[ ii ].push_back( 1.0 );
			packing_percentage_[ ii ].push_back( 0.0 );
			npd_percentage_[ ii ].push_back( 0.0 );
		}
	}
#endif
	//TR << "Finished computing state energies" << std::endl;

}

void MPIMultistateFitnessFunction::instruct_daemons_to_keep_last_entity()
{
	for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
		utility::send_integer_to_node( ii, keep_rotamer_assignment_for_last_entity );
	}
	if ( daemon_set() ) {
		daemon_set()->mark_last_entity_as_important();
	}
}

void MPIMultistateFitnessFunction::instruct_daemons_to_drop_entity( Entity const & entity )
{
	for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
		utility::send_integer_to_node( ii, discard_old_entity );
	}
	broadcast_entity_string( entity );
	if ( daemon_set() ) {
		daemon_set()->mark_entity_as_unimportant( entity );
	}
}

std::list< std::pair< MultistateFitnessFunction::Size, MultistateFitnessFunction::PoseOP > >
MPIMultistateFitnessFunction::recover_poses_from_states(
	Entity const & entity,
	utility::vector1< core::Size > const & which_states
)
{
	for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
		utility::send_integer_to_node( ii, geneate_pose_from_old_state );
	}
	broadcast_entity_string( entity );
	utility::vector1< int > which_states_ints( which_states );
	for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
		utility::send_integers_to_node( ii, which_states_ints );
	}
	std::list< std::pair< Size, PoseOP > > return_pose_list;
	if ( daemon_set() ) {
		std::list< std::pair< Size, PoseOP > > my_poses = daemon_set()->
			retrieve_relevant_poses_for_entity( entity, which_states );
		return_pose_list.splice( return_pose_list.end(), my_poses );
	}
	for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
		int n_poses_from_ii = utility::receive_integer_from_node( ii );
		//std::cout << "Ready to receive " << n_poses_from_ii << " pdb strings from node " << ii << std::endl;
		std::string pseudo_pdbname = "MPI_pdb_from_node_" + utility::to_string( ii );
		for ( int jj = 1; jj <= n_poses_from_ii; ++jj ) {
			//std::cout << "Receiving pose " << jj << " of " << n_poses_from_ii << std::flush;
			Size pack_daemon_index = utility::receive_integer_from_node( ii );
			std::string pdb_string = utility::receive_string_from_node( ii );

			PoseOP pose( new Pose );
			core::import_pose::pose_from_pdbstring( *pose, pdb_string, pseudo_pdbname );

			return_pose_list.emplace_back( pack_daemon_index, pose );
			//std::cout << "... done" << std::endl;
		}
	}
	return return_pose_list;
}

void MPIMultistateFitnessFunction::broadcast_entity_string( Entity const & entity )
{
	std::ostringstream oss;
	entity.write_checkpoint( oss );
	std::string entity_string = oss.str();

	for ( Size ii = 1; ii < MPI_nprocs_; ++ii ) {
		utility::send_string_to_node( ii, entity_string );
	}
}


}
}

