// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/calculators/ResidueCentralityCalculator.hh
/// @brief header file for ResidueCentralityCalculator class.
/// Roughly, fragment quality is number of fragments which are close to a pose in rmsd
/// @detailed
/// @author Tom Linsky (tlinsky@gmail.com)


#ifndef INCLUDED_devel_denovo_design_calculators_ResidueCentralityCalculator_hh
#define INCLUDED_devel_denovo_design_calculators_ResidueCentralityCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

//// C++ headers

// Parser headers
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

namespace devel {
namespace denovo_design {
namespace calculators {

class Node;
class Edge;
typedef utility::pointer::owning_ptr< Node > NodeOP;
typedef utility::pointer::owning_ptr< Edge > EdgeOP;

// Nodes and edges for dijkstra's algorithm
class Node : public utility::pointer::ReferenceCount
{
public:
	Node(std::string const id, core::Size const resi)
		: resi(resi), id(id), previous(NULL),
			distanceFromStart(9999),
			index( -1 ), lowlink( -1 )
	{
	}

public:
	core::Size resi;
	std::string id;
	NodeOP previous;
	int distanceFromStart;
	int index;
	int lowlink;
};

class Edge : public utility::pointer::ReferenceCount
{
public:
	Edge(NodeOP node1, NodeOP node2, int distance)
		: node1(node1), node2(node2), distance(distance)
	{
	}

	bool Connects(NodeOP node1, NodeOP node2)
	{
		return ( ( node1 == this->node1 && node2 == this->node2 ) ||
						 ( node1 == this->node2 && node2 == this->node1 ) );
	}

	bool contains( NodeOP node ) {
		return ( ( this->node1 == node ) || ( this->node2 == node ) );
	}

public:
	NodeOP node1;
	NodeOP node2;
	int distance;
};

NodeOP
ExtractSmallest(utility::vector1<NodeOP>& nodes);

utility::vector1<NodeOP>
AdjacentRemainingNodes(utility::vector1<NodeOP> const & nodes, utility::vector1<EdgeOP> const & edges, NodeOP node);

int
Distance(utility::vector1<EdgeOP> const & edges, NodeOP node1, NodeOP node2);

bool
Contains(utility::vector1<NodeOP> const & nodes, NodeOP node);

utility::vector1< utility::vector1<bool> >
neighbors( core::pose::Pose const & pose );


class ResidueCentralityCalculator : public core::pose::metrics::StructureDependentCalculator {
public:// constructor/destructor

	/// @brief default constructor
	ResidueCentralityCalculator();

	/// @brief copy constructor
	ResidueCentralityCalculator( ResidueCentralityCalculator const & rval );

	/// @brief destructor
	virtual ~ResidueCentralityCalculator();

public:// virtual constructor
	/// @brief make clone
  core::pose::metrics::PoseMetricCalculatorOP
	clone() const { return new ResidueCentralityCalculator( *this ); }

public:// mutators

protected:
  virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
  virtual std::string print( std::string const & key ) const;
  virtual void recompute( core::pose::Pose const & this_pose );

private:// private functions
	/// @brief determine the connectivity index of a residue
	core::Real
	connectivity_index( utility::vector1< NodeOP > const & nodes,
											utility::vector1< EdgeOP > const & edges,
											core::Size const resi ) const;

private:// member variables

	/// @brief the pose that recompute() was last called on
	core::pose::PoseOP last_pose_;

	/// @brief the computed list of centrality indices, by residue
	utility::vector1< core::Real > centralities_;

	/// @brief a list of residues excluded from the calculation
	utility::vector1< core::Size > excluded_residues_;

}; //ResidueCentralityCalculator


} // ns PoseMetricCalculators
} // ns denovo_design
} // ns devel

#endif
