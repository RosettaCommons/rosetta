// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pack_daemon/MultistateFitnessFunction.hh
/// @brief  declaration for class MultistateFitnessFunction to work with the PackDeamon classes
///         (not to be confused with J. Ashworth's MultiStateFitnessFunction class)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_pack_daemon_MultistateFitnessFunction_hh
#define INCLUDED_protocols_pack_daemon_MultistateFitnessFunction_hh

// Unit headers
#include <protocols/pack_daemon/MultistateFitnessFunction.fwd.hh>
#include <protocols/pack_daemon/MultistateAggregateFunction.fwd.hh>

// Package headers
#include <protocols/pack_daemon/PackDaemon.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/genetic_algorithm/FitnessFunction.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/vector0.hh>

// C++ headers
#include <list>

//Auto Headers
namespace protocols {
namespace pack_daemon {


class TopEntitySet : public utility::pointer::ReferenceCount
{
public:
	typedef protocols::genetic_algorithm::Entity          Entity;
	typedef protocols::genetic_algorithm::EntityOP        EntityOP;
	typedef utility::vector1< core::Real >                StateEnergies;
	typedef std::pair< StateEnergies, StateEnergies >     StateEnergiesAndNPDs;
	typedef std::pair< EntityOP, StateEnergiesAndNPDs >   EntityAndScore;
	typedef utility::vector1< EntityAndScore >            EntityHistory;
public:
	TopEntitySet();

	core::Size size() const;
	EntityAndScore const & operator[] ( core::Size index ) const;
	EntityAndScore & operator[] ( core::Size index );

	void desired_entity_history_size( core::Size setting );
	core::Size desired_entity_history_size() const;
	void clear();

	/// @brief Update the internal history after receiving a new entity;
	/// returns the list of entities that are to be discarded.  The boolean
	/// added_new_entity is set to true if the given entity was added to the
	/// the list of top entities.
	std::list< genetic_algorithm::EntityOP >
	update_entity_history(
		Entity const & ent,
		StateEnergiesAndNPDs const & seanpds,
		bool & added_new_entity
	);

	core::Size
	index_of_entity( Entity const & ent ) const;

	/// @brief remove the worst entity from the set and return it
	EntityAndScore pop();

private:
	Size          desired_entity_history_size_;
	Size          n_tied_for_worst_;
	EntityHistory top_entities_;
};

class MultistateFitnessFunction : public protocols::genetic_algorithm::FitnessFunction
{
public:
	typedef protocols::genetic_algorithm::FitnessFunction parent;
	typedef protocols::genetic_algorithm::Entity          Entity;
	typedef protocols::genetic_algorithm::EntityOP        EntityOP;
	typedef utility::vector1< core::Real >                StateEnergies;
	typedef std::pair< StateEnergies, StateEnergies >     StateEnergiesAndNPDs;
	typedef core::pose::PoseOP                            PoseOP;
	typedef core::pose::Pose                              Pose;
	typedef core::Real                                    Real;
	typedef core::Size                                    Size;

	typedef std::pair< EntityOP, StateEnergiesAndNPDs >   EntityAndScore;
	typedef utility::vector1< EntityAndScore >            EntityHistory;

public:
	MultistateFitnessFunction();
	virtual ~MultistateFitnessFunction();

	virtual core::Real evaluate( Entity & entity );

	StateEnergies const & state_energies() const;
	StateEnergies const & npd_properties() const;

	void daemon_set( DaemonSetOP ds );

	void aggregate_function( MultistateAggregateFunctionOP func );

	DaemonSetCOP daemon_set() const;
	MultistateAggregateFunctionCOP aggregate_function() const;

	void set_history_size( core::Size history_size );
	void clear_history();

	std::list< std::pair< Size, PoseOP > >
	recover_relevant_poses_for_entity( Entity const & );

protected:
	virtual void compute_state_energies( Entity const & entity );
	virtual core::Real compute_aggregate_score( Entity const & entity );

	virtual void instruct_daemons_to_keep_last_entity();
	virtual void instruct_daemons_to_drop_entity( Entity const & entity );

	virtual
	std::list< std::pair< Size, PoseOP > >
	recover_poses_from_states(
		Entity const &,
		utility::vector1< core::Size > const & which_states
	);

	// Read and write access for derived classes
	utility::vector1< core::Real > & state_energies();
	utility::vector1< core::Real > & npd_properties();

	DaemonSetOP daemon_set();
	MultistateAggregateFunctionOP aggregate_function();

	Size which_top_entity( Entity const & ent ) const;

	void update_entity_history( Entity const & ent );

private:

	DaemonSetOP daemon_set_;
	StateEnergies state_energies_;
	StateEnergies npd_properties_;
	MultistateAggregateFunctionOP aggregate_;

	TopEntitySet  top_entities_;

};

class MPIMultistateFitnessFunction : public MultistateFitnessFunction
{
public:
	MPIMultistateFitnessFunction();
	virtual ~MPIMultistateFitnessFunction();

	/// @brief Inform this MPIMultistateFitnessFunction how many PackDaemons
	/// are running on all nodes (counting this one) which it may be unaware of.
	void set_num_pack_daemons( Size n_daemons );

	/// @brief Inform this MPIMultistateFitnessFunction how many non-pairwise decomposable
	/// properties will be computed for all states
	void set_num_npd_properties( Size n_npd_properties );


	/// @brief Spin down the other nodes.  No entity evaluation may follow
	/// the spin down call
	void send_spin_down_signal();

#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	void print_load_balance_statistics( std::ostream & ) const;
	void reset_load_balance_statistics();
#endif

protected:

	/// @brief Override the base class implementation by broadcasting the
	/// entity to all slave nodes, and waiting for them to report the energies
	/// from their PackDaemons for that entity.  If this node has PackDaemons,
	/// then it will evaluate their energies, too.
	virtual void compute_state_energies( Entity const & entity );

	/// @brief Broadcast the instruction to DaemonSets on all other nodes to keep
	/// the entity that was just evaluated.
	virtual void instruct_daemons_to_keep_last_entity();

	/// @brief Broadcast the instruction to DaemonSets on all other nodes to discard
	/// an old entity which had been previously held on to.
	virtual void instruct_daemons_to_drop_entity( Entity const & entity );

	//// @brief Broadcast the entity and the set of states that correspond to
	/// a set of poses that are important.  The entity had better still live
	/// on the distant nodes.
	virtual
	std::list< std::pair< Size, PoseOP > >
	recover_poses_from_states(
		Entity const &,
		utility::vector1< core::Size > const & which_states
	);

private:

	void broadcast_entity_string( Entity const & entity );

private:
	Size MPI_nprocs_;

#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	utility::vector0< std::list< core::Real > > utilization_by_node_;
	utility::vector0< std::list< core::Real > > packing_percentage_;
	utility::vector0< std::list< core::Real > > npd_percentage_;
#endif

	//DaemonSetOP local_daemon_set_;

};


}
}

#endif
