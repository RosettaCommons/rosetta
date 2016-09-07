// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/TwelveANeighborGraph.hh
/// @brief  Twelve Angstrom Neighbor Graph class declaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_TwelveANeighborGraph_hh
#define INCLUDED_core_scoring_TwelveANeighborGraph_hh

// Unit Headers
#include <core/scoring/TwelveANeighborGraph.fwd.hh>

// Package Headers
#include <core/scoring/ContextGraph.hh>

// Project Headers
#include <core/graph/Graph.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {

class TwelveANeighborNode : public graph::Node
{
public:
	typedef graph::Node parent;

public:
	~TwelveANeighborNode() override;
	TwelveANeighborNode( graph::Graph*, Size node_id );

	void copy_from( Node const * source ) override;

	Size count_static_memory() const override;
	Size count_dynamic_memory() const override;
private:

};

class TwelveANeighborEdge : public graph::Edge
{
public:
	typedef graph::Edge parent;

public:
	~TwelveANeighborEdge() override;
	TwelveANeighborEdge(graph::Graph* owner, Size first_node_ind, Size second_node_ind);

	void copy_from( graph::Edge const * source ) override;

	Size count_static_memory() const override;
	Size count_dynamic_memory() const override;

private:

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	TwelveANeighborEdge();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class TwelveANeighborGraph : public ContextGraph
{
public:
	typedef ContextGraph parent;

public:

	~TwelveANeighborGraph() override;

	TwelveANeighborGraph();
	TwelveANeighborGraph( Size num_nodes );
	TwelveANeighborGraph( TwelveANeighborGraph const & source );
	TwelveANeighborGraph & operator = ( TwelveANeighborGraph const & source );

	
	Distance
	neighbor_cutoff() const override;

	
	void
	conditionally_add_edge(
		Size lower_node_id,
		Size upper_node_id,
		DistanceSquared dsq
	) override;

	
	ContextGraphOP
	clone() const override;

	
	void update_from_pose(
		pose::Pose const & pose
	) override;

	void delete_edge( graph::Edge * edge ) override;

protected:

	Size count_static_memory() const override;
	Size count_dynamic_memory() const override;

	graph::Node * create_new_node( Size node_index ) override;
	graph::Edge * create_new_edge( Size index1, Size index2) override;
	graph::Edge * create_new_edge( graph::Edge const * example_edge ) override;

private:
	static Distance const twelveA_;
	static DistanceSquared const twelveA_squared_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_TwelveANeighborGraph )
#endif // SERIALIZATION


#endif
