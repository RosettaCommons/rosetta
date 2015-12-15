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
/// @details
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


static THREAD_LOCAL basic::Tracer TR( "devel.denovo_design.calculators.ResidueCentralityCalculator" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace devel {
namespace denovo_design {
namespace calculators {

/// @brief default constructor
ResidueCentralityCalculator::ResidueCentralityCalculator() :
	core::pose::metrics::StructureDependentCalculator()
{
	centralities_.clear();
	excluded_residues_.clear();
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
	} else {
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

void Dijkstras( std::list< NodeOP > & nodes )
{
	while ( !nodes.empty() ) { //size() > 0 ) {
		NodeOP smallest = ExtractSmallest( nodes );
		TR.Debug << "smallest has distance=" << smallest->distanceFromStart << " and neighbors=" << smallest->neighbors.size() << std::endl;
		if ( smallest->distanceFromStart > 10000 ) {
			TR.Debug << "Dump of node table and shortest distances:";
			for ( std::list< NodeOP >::const_iterator it = nodes.begin(); it != nodes.end(); ++it ) {
				TR.Debug << (*it)->resi << " " << (*it)->distanceFromStart << std::endl;
			}
			utility_exit();
		}
		std::list< NodeOP > const & adjacentNodes( AdjacentRemainingNodes( nodes, smallest ) );
		TR << "Nodes adjacent to " << smallest->resi << ": ";
		for ( std::list< NodeOP >::const_iterator it = adjacentNodes.begin(); it != adjacentNodes.end(); ++it ) {
			TR.Debug << " " << (*it)->resi;
			int distance = smallest->distanceFromStart + 1;
			if ( distance < (*it)->distanceFromStart ) {
				(*it)->distanceFromStart = distance;
			}
		}
		TR.Debug << std::endl;
	}
}

// Find the node with the smallest distance,
// remove it, and return it.
NodeOP
ExtractSmallest( std::list< NodeOP > & nodes )
{
	if ( nodes.empty() /*size() == 0 */ ) return NULL;
	std::list< NodeOP >::iterator smallest( nodes.begin() );
	for ( std::list< NodeOP >::iterator current = ++(nodes.begin()); current != nodes.end(); ++current ) {
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
AdjacentRemainingNodes( std::list< NodeOP > const &, NodeOP node )
{
	assert( ! node->in_list );
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
	for ( std::list< NodeOP >::const_iterator it = nodes.begin(); it != nodes.end(); ++it ) {
		if ( node == *it ) {
			return true;
		}
	}
	return false;
}

void
find_neighbors( core::pose::Pose const & pose,
	std::list< NodeOP > const & nodes )
{
	for ( std::list< NodeOP >::const_iterator res_it_1 = nodes.begin(); res_it_1 != nodes.end(); ++res_it_1 ) {
		(*res_it_1)->neighbors.clear();
	}

	// Residues are in contact if >=1 atom from each residue is less far than this
	core::Real const distance_threshold( 5.0 );
	core::Real const max_possible_dist( 12.5 + distance_threshold );

	for ( std::list< NodeOP >::const_iterator res_it_1 = nodes.begin(); res_it_1 != nodes.end(); ++res_it_1 ) {
		core::conformation::Residue const & res_target( pose.residue( (*res_it_1)->resi ) );

		std::list< NodeOP >::const_iterator res_it_start = res_it_1;
		for ( std::list< NodeOP >::const_iterator res_it_2 = (++res_it_start); res_it_2 != nodes.end(); ++res_it_2 ) {
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

std::list< NodeOP >
generate_nodes( core::pose::Pose const & pose )
{
	std::list<NodeOP> nodes;

	// populate lists of nodes and edges
	for ( core::Size i=1; i<=pose.total_residue(); i++ ) {
		std::string resname = pose.residue( i ).name3() + boost::lexical_cast<std::string>( i );
		TR.Debug << "Adding node " << resname << std::endl;
		nodes.push_back( NodeOP( new Node( resname, i ) ) );
	} // for each residue

	TR << "nodes size=" << nodes.size() << ", computing edges..."; TR.flush();

	// find neighbors
	find_neighbors( pose, nodes );

	TR << "done." << std::endl;

	return nodes;
}

core::Real
ResidueCentralityCalculator::connectivity_index( std::list< NodeOP > const & nodes,
	core::Size const resi ) const {
	// reset nodes
	for ( std::list< NodeOP >::const_iterator it = nodes.begin(); it != nodes.end(); ++it ) {
		(*it)->in_list = true;
		if ( (*it)->resi == resi ) {
			(*it)->distanceFromStart = 0;
		} else {
			(*it)->distanceFromStart = INT_MAX;
		}
	}

	// run dijkstra's algorithm
	std::list< NodeOP > nodecopy(nodes);
	Dijkstras( nodecopy );

	core::Real running_sum = 0.0;
	for ( std::list< NodeOP >::const_iterator it = nodes.begin(); it != nodes.end(); ++it ) {
		core::Size shortest_path = (*it)->distanceFromStart;
		TR.Debug << "Shortest path from " << resi << " to " << (*it)->resi << " is " << shortest_path << std::endl;
		running_sum += shortest_path;
	} // for each residue
	return ( nodes.size() - 1 ) / running_sum;
}

/// @brief recompute interaction network
void
ResidueCentralityCalculator::recompute( core::pose::Pose const & pose )
{
	// convert to centroid and compare to saved pose; if it is the same, don't do anything.
	centralities_.clear();
	excluded_residues_.clear();

	// generate lists of nodes/edges
	std::list< NodeOP > const & nodes( generate_nodes( pose ) );
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		// if this residue is ala or gly, flag it by scoring it as -1
		if ( pose.residue( i ).name3() == "GLY" || pose.residue( i ).name3() == "ALA" ) {
			centralities_.push_back( -1 );
			excluded_residues_.push_back( i );
		} else {
			centralities_.push_back( connectivity_index( nodes, i ) );
		}
	}
	// compute mean for standardization
	core::Real mean( 0.0 );
	for ( core::Size i=1; i<=centralities_.size(); ++i ) {
		if ( centralities_[i] >= 0 ) {
			mean += centralities_[i];
		}
	}
	TR.Info << "Centralities size = " << centralities_.size() << "; excluded_residues size = " << excluded_residues_.size() << std::endl;
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


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
devel::denovo_design::calculators::ResidueCentralityCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( CEREAL_NVP( centralities_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( excluded_residues_ ) ); // utility::vector1<core::Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
devel::denovo_design::calculators::ResidueCentralityCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( centralities_ ); // utility::vector1<core::Real>
	arc( excluded_residues_ ); // utility::vector1<core::Size>
}

SAVE_AND_LOAD_SERIALIZABLE( devel::denovo_design::calculators::ResidueCentralityCalculator );
CEREAL_REGISTER_TYPE( devel::denovo_design::calculators::ResidueCentralityCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( devel_denovo_design_calculators_ResidueCentralityCalculator )
#endif // SERIALIZATION
