// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/MultiplexedAnnealableGraph.cc
/// @brief  Multiplexing packing graph container.
/// @author Alex Ford (fordas@uw.edu)

#include <boost/foreach.hpp>

// Unit Headers
#include <core/pack/interaction_graph/MultiplexedAnnealableGraph.hh>

/// Utility headers
#include <utility/exit.hh>

//STL Headers
#include <list>
#include <cassert>
#include <stdexcept>

//ObjexxFCL Headers
#include <ObjexxFCL/FArray1A.hh>

using namespace ObjexxFCL;

namespace core {
namespace pack {
namespace interaction_graph {

/// @brief Default constructor.
///
MultiplexedAnnealableGraph::MultiplexedAnnealableGraph():
	AnnealableGraphBase()
{}

/// @brief Copy constructor.
///
MultiplexedAnnealableGraph::MultiplexedAnnealableGraph(MultiplexedAnnealableGraph const & other) :
	AnnealableGraphBase(other),
	subgraphs(other.subgraphs)
{}

/// @brief Setup constructor.
///
MultiplexedAnnealableGraph::MultiplexedAnnealableGraph(SubgraphContainer const & target_subgraphs) :
	AnnealableGraphBase(),
	subgraphs(target_subgraphs)
{}

/// @brief Destructor.
///
MultiplexedAnnealableGraph::~MultiplexedAnnealableGraph()
{}


int MultiplexedAnnealableGraph::get_num_nodes() const
{
	if ( subgraphs.empty() ) {
		throw std::logic_error("MultiplexedAnnealableGraph with empty subgraphs." );
	}

	int num_nodes = subgraphs.front()->get_num_nodes();

	BOOST_FOREACH ( AnnealableGraphBaseOP subgraph, subgraphs ) {
		if ( subgraph->get_num_nodes() != num_nodes ) {
			throw std::logic_error("inconsistent num_nodes in subgraphs.");
		}
	}

	return num_nodes;
}

int MultiplexedAnnealableGraph::get_num_states_for_node(int node) const
{
	if ( subgraphs.empty() ) {
		throw std::logic_error("MultiplexedAnnealableGraph with empty subgraphs." );
	}

	int num_states = subgraphs.front()->get_num_states_for_node(node);

	BOOST_FOREACH ( AnnealableGraphBaseOP subgraph, subgraphs ) {
		if ( subgraph->get_num_states_for_node(node) != num_states ) {
			throw std::logic_error("inconsistent num_states_for_node in subgraphs.");
		}
	}

	return num_states;
}

int MultiplexedAnnealableGraph::get_num_total_states() const
{
	if ( subgraphs.empty() ) {
		throw std::logic_error("MultiplexedAnnealableGraph with empty subgraphs." );
	}

	int num_states = subgraphs.front()->get_num_total_states();

	BOOST_FOREACH ( AnnealableGraphBaseOP subgraph, subgraphs ) {
		if ( subgraph->get_num_total_states() != num_states ) {
			throw std::logic_error("inconsistent num_total_states in subgraphs.");
		}
	}

	return num_states;
}

void MultiplexedAnnealableGraph::prepare_for_simulated_annealing()
{
	BOOST_FOREACH ( AnnealableGraphBaseOP subgraph, subgraphs ) {
		subgraph->prepare_for_simulated_annealing();
	}
}

void MultiplexedAnnealableGraph::blanket_assign_state_0()
{
	BOOST_FOREACH ( AnnealableGraphBaseOP subgraph, subgraphs ) {
		subgraph->blanket_assign_state_0();
	}
}

bool MultiplexedAnnealableGraph::any_vertex_state_unassigned() const
{
	BOOST_FOREACH ( AnnealableGraphBaseOP subgraph, subgraphs ) {
		if ( subgraph->any_vertex_state_unassigned() ) {
			return true;
		}
	}

	return false;
}

core::PackerEnergy MultiplexedAnnealableGraph::set_state_for_node(int node_ind, int new_state)
{
	core::PackerEnergy result = 0;

	BOOST_FOREACH ( AnnealableGraphBaseOP subgraph, subgraphs ) {
		result += subgraph->set_state_for_node(node_ind, new_state);
	}

	return result;
}

core::PackerEnergy MultiplexedAnnealableGraph::set_network_state( ObjexxFCL::FArray1_int & node_states )
{
	core::PackerEnergy result = 0;

	BOOST_FOREACH ( AnnealableGraphBaseOP subgraph, subgraphs ) {
		result += subgraph->set_network_state(node_states);
	}

	return result;
}

void MultiplexedAnnealableGraph::consider_substitution(
	int node_ind,
	int new_state,
	core::PackerEnergy & delta_energy,
	core::PackerEnergy & prev_energy_for_node)
{
	delta_energy = 0;
	prev_energy_for_node = 0;

	BOOST_FOREACH ( AnnealableGraphBaseOP subgraph, subgraphs ) {
		core::PackerEnergy subgraph_delta = 0;
		core::PackerEnergy subgraph_prev = 0;
		subgraph->consider_substitution( node_ind, new_state, subgraph_delta, subgraph_prev);

		delta_energy += subgraph_delta;
		prev_energy_for_node += subgraph_prev;
	}
}

core::PackerEnergy MultiplexedAnnealableGraph::commit_considered_substitution()
{
	core::PackerEnergy result = 0;

	BOOST_FOREACH ( AnnealableGraphBaseOP subgraph, subgraphs ) {
		result += subgraph->commit_considered_substitution();
	}

	return result;
}

core::PackerEnergy MultiplexedAnnealableGraph::get_energy_current_state_assignment()
{
	core::PackerEnergy result = 0;

	BOOST_FOREACH ( AnnealableGraphBaseOP subgraph, subgraphs ) {
		result += subgraph->get_energy_current_state_assignment();
	}

	return result;
}

void MultiplexedAnnealableGraph::set_errorfull_deltaE_threshold(core::PackerEnergy deltaE)
{
	BOOST_FOREACH ( AnnealableGraphBaseOP subgraph, subgraphs ) {
		subgraph->set_errorfull_deltaE_threshold(deltaE);
	}
}

}
}
}

