// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ContextGraph.hh
/// @brief  Context graph class declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_ContextGraph_hh
#define INCLUDED_core_scoring_ContextGraph_hh

// Unit Headers
#include <core/scoring/ContextGraph.fwd.hh>

// Project Headers
#include <core/graph/Graph.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

class ContextGraph : public graph::Graph {
public:
	typedef graph::Graph parent;

public:

	virtual ~ContextGraph();

	ContextGraph();
	ContextGraph(Size num_nodes);
	ContextGraph( ContextGraph const & source );
	ContextGraph & operator = ( ContextGraph const & source );

	virtual
	Distance
	neighbor_cutoff() const = 0;


	virtual
	void
	conditionally_add_edge(
		Size lower_node_id,
		Size upper_node_id,
		DistanceSquared dsq
	) = 0;

	virtual
	ContextGraphOP
	clone() const = 0;

	virtual
	void update_from_pose(
		pose::Pose const & pose
	) = 0;

protected:
	virtual Size count_static_memory() const = 0;
	virtual Size count_dynamic_memory() const;


};

} // scoring
} // core

#endif
