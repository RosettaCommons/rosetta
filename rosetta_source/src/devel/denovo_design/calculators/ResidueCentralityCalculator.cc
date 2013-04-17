// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/denovo_design/calculators/ResidueCentralityCalculator.cc
/// @brief calculate centrality for all (or a subset of) residues
/// @detailed
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit Headers
#include <devel/denovo_design/calculators/ResidueCentralityCalculator.hh>

// Project Headers
#include <basic/MetricValue.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <boost/lexical_cast.hpp>

//// C++ headers
#include <ObjexxFCL/format.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>


static basic::Tracer TR("devel.denovo_design.calculators.ResidueCentralityCalculator");

namespace devel {
namespace denovo_design {
namespace calculators {

/// @brief default constructor
ResidueCentralityCalculator::ResidueCentralityCalculator() :
	core::pose::metrics::StructureDependentCalculator(),
	last_pose_( NULL )
{
	centralities_.clear();
	excluded_residues_.clear();
}

/// @brief copy constructor
ResidueCentralityCalculator::ResidueCentralityCalculator( ResidueCentralityCalculator const & rval ) :
	core::pose::metrics::StructureDependentCalculator( rval ),
	last_pose_( NULL ),
	centralities_( rval.centralities_ ),
	excluded_residues_( rval.excluded_residues_ )
{
	if ( rval.last_pose_ ) {
		last_pose_ = new core::pose::Pose( *(rval.last_pose_) );
	}
}

/// @brief destructor
ResidueCentralityCalculator::~ResidueCentralityCalculator(){}


/// @brief
void
ResidueCentralityCalculator::lookup( std::string const & key,
																		 basic::MetricValueBase * valptr ) const
{
	if ( key == "centrality" ) {
		basic::check_cast( valptr, &centralities_, "residue centrality list" );
		(static_cast< basic::MetricValue<utility::vector1< core::Real > > *>(valptr))->set( centralities_ );
	} else if ( key == "excluded_residues" ) {
		basic::check_cast( valptr, &excluded_residues_, "list of excluded residues" );
		(static_cast< basic::MetricValue<utility::vector1< core::Size > > *>(valptr))->set( excluded_residues_ );
	}	else {
		TR << "ResidueCentralityCalculator cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}

} //lookup


// @brief
std::string
ResidueCentralityCalculator::print( std::string const & key ) const
{
	std::string result;
	if ( key == "centrality" ) {
		result += "[ ";
		for ( core::Size i=1; i<=centralities_.size(); ++i ) {
			result += boost::lexical_cast<std::string>( centralities_[i] ) + " ";
		}
		result += "]";
  } else if ( key == "excluded_residues" ) {
		result += "[ ";
		for ( core::Size i=1; i<=excluded_residues_.size(); ++i ) {
			result += boost::lexical_cast<std::string>( excluded_residues_[i] ) + " ";
		}
		result += "]";
	} else {
		basic::Error() << "ResidueCentralityCalculator cannot compute metric " << key << std::endl;
	}
	return result;
} // apply

void Dijkstras(utility::vector1<NodeOP> &nodes, utility::vector1<EdgeOP> const & edges)
{
	while (nodes.size() > 0) {
		NodeOP smallest = ExtractSmallest(nodes);
		utility::vector1<NodeOP> const & adjacentNodes( AdjacentRemainingNodes(nodes, edges, smallest) );
		core::Size const size = adjacentNodes.size();
		for (core::Size i=1; i<=size; ++i) {
			int distance = Distance(edges, smallest, adjacentNodes[i] ) +
				smallest->distanceFromStart;
			if (distance < adjacentNodes[i]->distanceFromStart)	{
				adjacentNodes[i]->distanceFromStart = distance;
				adjacentNodes[i]->previous = smallest;
			}
		}
	}
}

// Find the node with the smallest distance,
// remove it, and return it.
NodeOP
ExtractSmallest(utility::vector1<NodeOP>& nodes)
{
	core::Size const size = nodes.size();
	if (size == 0) return NULL;
	core::Size smallestPosition = 1;
	NodeOP smallest = nodes.at(1);
	for (core::Size i=2; i<=size; ++i) {
		NodeOP current = nodes.at(i);
		if (current->distanceFromStart <
				smallest->distanceFromStart) {
			smallest = current;
			smallestPosition = i;
		}
	}
	nodes.erase(nodes.begin() + smallestPosition - 1);
	return smallest;
}

// Return all nodes adjacent to 'node' which are still
// in the 'nodes' collection.
utility::vector1<NodeOP>
AdjacentRemainingNodes(utility::vector1<NodeOP> const & nodes, utility::vector1<EdgeOP> const & edges, NodeOP node)
{
	utility::vector1<NodeOP> adjacentNodes;
	core::Size const size = edges.size();
	for(core::Size i=1; i<=size; ++i) {
		EdgeOP edge = edges.at(i);
		NodeOP adjacent = NULL;
		if (edge->node1 == node) {
			adjacent = edge->node2;
		}
		else if (edge->node2 == node) {
			adjacent = edge->node1;
		}
		if (adjacent && Contains(nodes, adjacent)) {
			adjacentNodes.push_back(adjacent);
		}
	}
	return adjacentNodes;
}

// Return distance between two connected nodes
int
Distance(utility::vector1<EdgeOP> const & edges, NodeOP node1, NodeOP node2)
{
	core::Size const size = edges.size();
	for(core::Size i=1; i<=size; ++i)
		{
			EdgeOP edge = edges.at(i);
			if (edge->Connects(node1, node2))
				{
					return edge->distance;
				}
		}
	return -1; // should never happen
}

// Does the 'nodes' vector contain 'node'
bool
Contains(utility::vector1<NodeOP> const & nodes, NodeOP node)
{
	core::Size const size = nodes.size();
	for(core::Size i=1; i<=size; ++i)
		{
			if (node == nodes.at(i))
				{
					return true;
				}
		}
	return false;
}

utility::vector1< utility::vector1< bool > >
neighbors( core::pose::Pose const & pose )
{
	// Residues are in contact if >=1 atom from each residue is less far than this
	core::Real const distance_threshold( 5.0 );
	core::Real const max_possible_dist( 12.5 + distance_threshold );

	utility::vector1<bool> vecn( pose.total_residue(), false );
	utility::vector1< utility::vector1< bool > > vec( pose.total_residue(), vecn );

	for ( core::Size resi=1; resi<=pose.total_residue()-1; ++resi ) {
		core::conformation::Residue const & res_target( pose.residue( resi ) );

		core::Size resi_neighbors( 0 );
		for ( core::Size i=resi+1; i<=pose.total_residue(); ++i ) {
			core::conformation::Residue const & resj( pose.residue( i ) );
			core::Real closest( -1 );

			// if the distance is less than the maximum possible interaction distance (arg side chain + arg side chain ~= 12.5 A)
			// then we compute distances for each atom
			if ( res_target.xyz( res_target.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) <= max_possible_dist ) {
				core::conformation::Atoms const & atoms1( res_target.atoms() );
				core::conformation::Atoms const & atoms2( resj.atoms() );
				for ( core::Size a1=1; a1<=atoms1.size(); ++a1 ) {
					for ( core::Size a2=1; a2<=atoms2.size(); ++a2 ) {
						core::Real const distance( atoms1[a1].xyz().distance( atoms2[a2].xyz() ) );
						TR.Debug << "r1=" << resi << ", a1=" << atoms1[a1] << ", r2=" << i << ", a2=" << atoms2[a2] << ", dist=" << distance << std::endl;
						if ( closest <= 0 || closest > distance ) {
							closest = distance;
						}
						// stop if any atom is within the distance threshold
						if ( closest <= distance_threshold ) {
							break;
						}
					}
					if ( closest <= distance_threshold ) {
						break;
					}
				}
				if( closest <= distance_threshold ){
					vec[resi][i] = true;
					vec[i][resi] = true;
					++resi_neighbors;
					TR.Debug<<"Adding residue "<<i<<" as a neighbor of "<<resi<<" -- Count = "<<resi_neighbors<<std::endl;
				}
			}
		}
	}
	return vec;
}

utility::vector1<NodeOP>
generate_nodes( core::pose::Pose const & pose )
{
	utility::vector1<NodeOP> nodes;
	// populate lists of nodes and edges
	core::Size nligands = 0;
	for ( core::Size i=1; i<=pose.total_residue(); i++ ){
		std::string resname = pose.residue( i ).name3() + boost::lexical_cast<std::string>( i );
		TR.Debug << "Adding node " << resname << std::endl;
		NodeOP newnode = new Node( resname, i );
		nodes.push_back( newnode );
	} // for each residue
	TR << "nodes size=" << nodes.size() << ", computing edges..." << std::endl;
	return nodes;
}

utility::vector1<EdgeOP>
generate_edges( core::pose::Pose const & pose, utility::vector1<NodeOP> const & nodes )
{
	runtime_assert( pose.total_residue() == nodes.size() );
	utility::vector1<EdgeOP> edges;

	// get neighbor list
	utility::vector1< utility::vector1< bool > > const & n_vec( neighbors( pose ) );
	runtime_assert( n_vec.size() == pose.total_residue() );
	runtime_assert( n_vec[1].size() == pose.total_residue() );

	// populate lists of nodes and edges
	core::Size nligands = 0;
	for ( core::Size i=1; i<=pose.total_residue()-1; ++i ){
		NodeOP newnode = nodes[i];
		std::string resname = pose.residue( i ).name3() + boost::lexical_cast<std::string>( i );

		for ( core::Size j=i+1; j<=pose.total_residue(); ++j ){
			if ( n_vec[i][j] ){
				TR.Debug << "Adding edge from " << nodes[j]->id << " to " << newnode->id << std::endl;
				edges.push_back( new Edge( nodes[j], newnode, 1 ) );
			}
		} // for neighbors
	} // for each residue
	TR << "edges size=" << edges.size() << " and total residues in pose=" << pose.total_residue() << std::endl;
	return edges;
}

core::Real
ResidueCentralityCalculator::connectivity_index( utility::vector1< NodeOP> const & nodes,
																								 utility::vector1< EdgeOP> const & edges,
																								 core::Size const resi ) const {
	// reset nodes
	for ( core::Size i=1; i<= nodes.size(); ++i ) {
		if ( i == resi ) {
			nodes[i]->distanceFromStart = 0;
		} else {
			nodes[i]->distanceFromStart = INT_MAX;
		}

	}
	// run dijkstra's algorithm
	utility::vector1<NodeOP> nodecopy(nodes);
	Dijkstras( nodecopy, edges );

	core::Real running_sum = 0.0;
	for ( core::Size i=1; i<=nodes.size(); i++ ) {
		core::Size shortest_path = nodes[i]->distanceFromStart;
		TR.Debug << "Shortest path from " << resi << " to " << i << " is " << shortest_path << std::endl;
		running_sum += shortest_path;
	} // for each residue
	return ( nodes.size() - 1 ) / running_sum;
}

/// @brief recompute interaction network
void
ResidueCentralityCalculator::recompute( core::pose::Pose const & pose )
{
	// convert to centroid and compare to saved pose; if it is the same, don't do anything.

	last_pose_ = new core::pose::Pose( pose );
	centralities_.clear();
	// generate lists of nodes/edges
	utility::vector1<NodeOP> const & nodes( generate_nodes( pose ) );
	utility::vector1<EdgeOP> const & edges( generate_edges( pose , nodes ) );
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		// if this residue is ala or gly, flag it by scoring it as -1
		if ( pose.residue( i ).name3() == "GLY" || pose.residue( i ).name3() == "ALA" ) {
			centralities_.push_back( -1 );
			excluded_residues_.push_back( i );
		} else {
			centralities_.push_back( connectivity_index( nodes, edges, i ) );
		}
	}
	// compute mean for standardization
	core::Real mean( 0.0 );
	for ( core::Size i=1; i<=centralities_.size(); ++i ) {
		if ( centralities_[i] >= 0 ) {
			mean += centralities_[i];
		}
	}
	mean /= ( centralities_.size() - excluded_residues_.size() );

	// compute std dev
	core::Real dev( 0.0 );
	for ( core::Size i=1; i<=centralities_.size(); ++i ) {
		if ( centralities_[i] >= 0 ) {
			dev += ( centralities_[i] - mean ) * ( centralities_[i] - mean );
		}
	}
	dev /= ( centralities_.size() - 1 - excluded_residues_.size() );
	dev = std::sqrt( dev );

	TR << "mean=" << mean << ", dev=" << dev << std::endl;

	// compute z-scores
	for ( core::Size i=1; i<=centralities_.size(); ++i ) {
		TR << "Centrality for residue " << i << " : " << centralities_[i] << " : ";
		if ( centralities_[i] >= 0 ) {
			centralities_[i] = ( centralities_[i] - mean ) / dev;
		}
		TR << centralities_[i] << std::endl;
	}
}

} // calculators
} // denovo_design
} // devel

