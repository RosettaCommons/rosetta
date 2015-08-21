// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/TenANeighborGraph.hh
/// @brief  Ten Angstrom Neighbor Graph class declaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_TenANeighborGraph_hh
#define INCLUDED_core_scoring_TenANeighborGraph_hh

// Unit Headers
#include <core/scoring/TenANeighborGraph.fwd.hh>

// Package Headers
#include <core/scoring/ContextGraph.hh>

// Project Headers
#include <core/graph/Graph.hh>

#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

class TenANeighborNode : public graph::Node
{
public:
	typedef graph::Node parent;

public:
	virtual ~TenANeighborNode();
	TenANeighborNode( graph::Graph*, Size node_id );

	virtual void copy_from( Node const * source );

	virtual Size count_static_memory() const;
	virtual Size count_dynamic_memory() const;

	/// @brief set the neighbor mass for a vertex
	void neighbor_mass( Real mass );

	/// @brief return the neighbor mass for a vertex
	inline
	Real neighbor_mass() const
	{ return neighbor_mass_; }

	inline
	Real sum_of_neighbors_masses() const
	{
		if ( since_last_sonm_update_ == TIME_TO_UPDATE ) {
			update_neighbor_mass_sum();
		}
		return sum_of_neighbors_masses_;
	}

	/// @brief To be called by TenANeighborEdge only
	/// called at the time of edge addition
	void add_neighbors_mass( Real neighbors_mass ) {
		sum_of_neighbors_masses_ += neighbors_mass;
		++since_last_sonm_update_;
	}

	/// @brief To be called by TenANeighborEdge only
	/// called at the time of edge deletion
	void subtract_neighbors_mass( Real neighbors_mass ) {
		sum_of_neighbors_masses_ -= neighbors_mass;
		++since_last_sonm_update_;
	}
private:
	void
	update_neighbor_mass_sum() const;

private:
	static Size const TIME_TO_UPDATE = 1024;

	Real neighbor_mass_;
	mutable Real sum_of_neighbors_masses_;
	mutable Size since_last_sonm_update_;
};

class TenANeighborEdge : public graph::Edge
{
public:
	typedef graph::Edge parent;

public:
	virtual ~TenANeighborEdge();
	TenANeighborEdge(graph::Graph* owner, Size first_node_ind, Size second_node_ind);

	virtual void copy_from( graph::Edge const * source );

	virtual Size count_static_memory() const;
	virtual Size count_dynamic_memory() const;

private:
	inline
	TenANeighborNode const *
	get_TenANode( Size node_index ) const {
		debug_assert( dynamic_cast< TenANeighborNode const * > ( get_node( node_index ) ) );
		return static_cast< TenANeighborNode const * > ( get_node( node_index ) );
	}

	inline
	TenANeighborNode *
	get_TenANode( Size node_index ) {
		debug_assert( dynamic_cast< TenANeighborNode * > ( get_node( node_index ) ) );
		return static_cast< TenANeighborNode * > ( get_node( node_index ) );
	}

};

class TenANeighborGraph : public ContextGraph
{
public:
	typedef ContextGraph parent;

public:

	virtual ~TenANeighborGraph();

	TenANeighborGraph();
	TenANeighborGraph(Size num_nodes);
	TenANeighborGraph( TenANeighborGraph const & source );
	TenANeighborGraph & operator = ( TenANeighborGraph const & source );

	virtual
	Distance
	neighbor_cutoff() const;

	virtual
	void
	conditionally_add_edge(
		Size lower_node_id,
		Size upper_node_id,
		DistanceSquared dsq
	);

	virtual
	ContextGraphOP
	clone() const;

	virtual
	void update_from_pose(
		pose::Pose const & pose
	);

	virtual void delete_edge( graph::Edge * edge );


protected:

	virtual Size count_static_memory() const;
	virtual Size count_dynamic_memory() const;

	virtual graph::Node * create_new_node( Size node_index );
	virtual graph::Edge * create_new_edge( Size index1, Size index2);
	virtual graph::Edge * create_new_edge( graph::Edge const * example_edge );

private:
	static Distance const tenA_;
	static DistanceSquared const tenA_squared_;

	boost::unordered_object_pool< TenANeighborEdge > * tenA_edge_pool_;

};

} // scoring
} // core

#endif
