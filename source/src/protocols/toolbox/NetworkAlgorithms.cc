// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/NetworkAlgorithms.cc
/// @brief calculate various network properties from poses
/// @details
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit Headers
#include <protocols/toolbox/NetworkAlgorithms.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>

// Utility headers

//// C++ headers
#include <limits.h> // For INT_MAX

//Auto Headers


static basic::Tracer TR( "protocols.toolbox.NetworkAlgorithms" );

namespace protocols {
namespace toolbox {

/// @brief default constructor
ResidueNetwork::ResidueNetwork() = default;

/// @brief destructor
ResidueNetwork::~ResidueNetwork() = default;

core::Real
ResidueNetwork::connectivity_index( core::Size const resi ) const
{
	runtime_assert( nodes_.size() > 0 );

	// run dijkstra's algorithm
	dijkstras( resi );

	core::Real running_sum = 0.0;
	for ( auto const & node : nodes_ ) {
		core::Size shortest_path = node->distanceFromStart;
		TR.Debug << "Shortest path from " << resi << " to " << node->resi << " is " << shortest_path << std::endl;
		running_sum += shortest_path;
	} // for each residue
	return ( nodes_.size() - 1 ) / running_sum;
}

/// @brief calculates the average shortest path length of the network
core::Real
ResidueNetwork::average_shortest_path_length() const
{
	runtime_assert( nodes_.size() > 0 );

	core::Real total_path_length = 0.0;

	//iterate over all starting notes
	for ( auto const & it : nodes_ ) {
		dijkstras (it->resi);

		//add the paths from resi it to all other residues
		for ( auto const & node : nodes_ ) {
			core::Size shortest_path = node->distanceFromStart;
			total_path_length += shortest_path;
		}
	}

	total_path_length = ((total_path_length / nodes_.size()) / (nodes_.size()-1));
	TR << "AveragePathLength " << total_path_length << std::endl;
	return total_path_length;
}

/// @brief run Dijkstra's shortest path algorithm on the given list of nodes
/// after execution, the "distanceFromStart" variable of each node will contain the distance from residue resi
void
ResidueNetwork::dijkstras( core::Size const resi ) const
{
	// reset nodes
	for ( auto const & node : nodes_ ) {
		node->in_list = true;
		if ( node->resi == resi ) {
			node->distanceFromStart = 0;
		} else {
			node->distanceFromStart = INT_MAX;
		}
	}

	// copy the node list so that the original doesn't get modified
	std::list< NodeOP > nodes( nodes_ );
	while ( !nodes.empty() ) {//size() > 0 ) {
		NodeOP smallest = ExtractSmallest( nodes );
		TR.Debug << "smallest has distance=" << smallest->distanceFromStart << " and neighbors=" << smallest->neighbors.size() << std::endl;
		if ( smallest->distanceFromStart > 10000 ) {
			TR.Debug << "Dump of node table and shortest distances:";
			for ( std::list< NodeOP >::const_iterator it = nodes.begin(); it != nodes.end(); ++it ) {
				TR.Debug << (*it)->resi << " " << (*it)->distanceFromStart << std::endl;
			}
			utility_exit_with_message("In ResidueNetwork, the smallest distance to start is too big.");
		}
		std::list< NodeOP > const & adjacentNodes( AdjacentRemainingNodes( smallest ) );
		TR.Debug << "Nodes adjacent to " << smallest->resi << ": ";
		for ( auto const & adjacentNode : adjacentNodes ) {
			TR.Debug << " " << adjacentNode->resi;
			int distance = smallest->distanceFromStart + 1;
			if ( distance < adjacentNode->distanceFromStart ) {
				adjacentNode->distanceFromStart = distance;
			}
		}
		TR.Debug << std::endl;
	}
}

/// @brief create a network from an input pose
void
ResidueNetwork::create_from_pose( core::pose::Pose const & pose )
{
	nodes_.clear();

	// populate lists of nodes and edges
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		std::string resname = pose.residue( i ).name3() + std::to_string( i );
		TR.Debug << "Adding node " << resname << std::endl;
		nodes_.push_back( utility::pointer::make_shared< Node >( resname, i ) );
	} // for each residue

	TR << "nodes size=" << nodes_.size() << ", computing edges..."; TR.flush();

	// find neighbors
	generate_edges( pose );

	TR << "done." << std::endl;
}

// Find the node with the smallest distance,
// remove it, and return it.
NodeOP
ExtractSmallest( std::list< NodeOP > & nodes )
{
	// AMW: cppcheck thinks this conditional might be inefficient
	if ( nodes.empty() /*size() == 0 */ ) return nullptr;
	auto smallest( nodes.begin() );
	for ( auto current = ++(nodes.begin()); current != nodes.end(); ++current ) {
		if ( (*current)->distanceFromStart < (*smallest)->distanceFromStart ) {
			smallest = current;
		}
	}
	NodeOP retval( *smallest );
	nodes.erase( smallest );
	retval->in_list = false;
	return retval;
}

// Return all nodes adjacent to 'node' which are still
// in the 'nodes' collection.
std::list< NodeOP >
AdjacentRemainingNodes( NodeOP node )
{
	debug_assert( ! node->in_list );
	std::list< NodeOP > adjacentNodes;
	for ( std::list< NodeOP >::const_iterator cur_neighbor = node->neighbors.begin();
			cur_neighbor != node->neighbors.end();
			++cur_neighbor ) {
		if ( (*cur_neighbor)->in_list ) {
			//if ( Contains( nodes, *cur_neighbor ) ) {
			adjacentNodes.push_back( *cur_neighbor );
		}
	}
	return adjacentNodes;
}

// Does the 'nodes' vector contain 'node'
bool
Contains( std::list< NodeOP > const & nodes, NodeCOP node)
{
	for ( auto const & it : nodes ) {
		if ( node == it ) {
			return true;
		}
	}
	return false;
}

/// @brief empties edges
void
ResidueNetwork::clear_edges()
{
	for ( auto const & res_it_1 : nodes() ) {
		res_it_1->neighbors.clear();
	}
}

void
DistanceResidueNetwork::generate_edges( core::pose::Pose const & pose )
{
	clear_edges();

	// Residues are in contact if >=1 atom from each residue is less far than this
	core::Real const distance_threshold( 5.0 );
	core::Real const max_possible_dist( 12.5 + distance_threshold );

	for ( auto res_it_1 = nodes().begin(); res_it_1 != nodes().end(); ++res_it_1 ) {
		core::conformation::Residue const & res_target( pose.residue( (*res_it_1)->resi ) );

		auto res_it_start = res_it_1;
		for ( auto res_it_2 = (++res_it_start); res_it_2 != nodes().end(); ++res_it_2 ) {
			core::conformation::Residue const & resj( pose.residue( (*res_it_2)->resi ) );
			core::Real closest( -1.0 );
			// if the distance is less than the maximum possible interaction distance (arg side chain + arg side chain ~= 12.5 A)
			// then we compute distances for each atom
			core::Real nbr_dist( res_target.xyz( res_target.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) );
			if ( nbr_dist <= distance_threshold ) {
				(*res_it_1)->neighbors.push_back( *res_it_2 );
				(*res_it_2)->neighbors.push_back( *res_it_1 );
				TR.Debug<<"Adding residue "<<(*res_it_2)->resi<<" as a neighbor of "<<(*res_it_1)->resi<<" by nbr. -- Count = "<<(*res_it_1)->neighbors.size()<<","<<(*res_it_2)->neighbors.size()<<std::endl;
			} else if ( nbr_dist <= max_possible_dist ) {
				core::conformation::Atoms const & atoms1( res_target.atoms() );
				core::conformation::Atoms const & atoms2( resj.atoms() );
				for ( core::Size a1=1; a1<=atoms1.size(); ++a1 ) {
					for ( core::Size a2=1; a2<=atoms2.size(); ++a2 ) {
						core::Real const distance( atoms1[a1].xyz().distance( atoms2[a2].xyz() ) );
						TR.Debug << "r1=" << (*res_it_1)->resi << ", a1=" << atoms1[a1] << ", r2=" << (*res_it_2)->resi << ", a2=" << atoms2[a2] << ", dist=" << distance << std::endl;
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
				if ( closest <= distance_threshold ) {
					(*res_it_1)->neighbors.push_back( *res_it_2 );
					(*res_it_2)->neighbors.push_back( *res_it_1 );
					TR.Debug<<"Adding residue "<<(*res_it_2)->resi<<" as a neighbor of "<<(*res_it_1)->resi<<" -- Count = "<<(*res_it_1)->neighbors.size()<<","<<(*res_it_2)->neighbors.size()<<std::endl;
				}
			}
		}
	}
}

void
CovalentResidueNetwork::generate_edges( core::pose::Pose const & pose )
{
	clear_edges();

	//identify covalent bonds between residues
	auto res_it_1 = nodes().begin();
	for ( core::Size i=1; i != pose.size(); ++i ) {
		auto res_it_2 = nodes().begin();
		std::advance(res_it_2,i);
		for ( core::Size j=i + 1; j != pose.size() + 1; ++j ) {
			if ( pose.residue(i).is_bonded(pose.residue(j)) ) {
				//TR << "bond found " << i << " " << j << std::endl;
				if ( j > i + 1 ) {
					TR << "Non-peptide bond found (probably disulfide) " << i << "-" << j << std::endl;
				}
				(*res_it_1)->neighbors.push_back(*(res_it_2));
				(*res_it_2)->neighbors.push_back(*(res_it_1));
			}
			++res_it_2;
		}
		++res_it_1;
	}
}


} // toolbox
} // protocols

