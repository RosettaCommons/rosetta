// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/PackerNeighborGraphFilter.cc
/// @brief
/// @details
///   Contains currently:
///
///
/// @author Florian Richter (floric@u.washington.edu ), march 2009

// Unit Headers
#include <protocols/simple_filters/PackerNeighborGraphFilter.hh>

// Package Headers

// Project Headers
#include <utility/graph/Graph.hh>
#include <core/pack/packer_neighbors.hh>
//#include <core/pose/Pose.hh>
#include <core/types.hh>


// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


//// C++ headers
static THREAD_LOCAL basic::Tracer tr( "protocols.filters.PackerNeighborGraphFilter" );

namespace protocols {
namespace simple_filters {


void
RegionalConnections::check_if_connected_residues_belong_to_regions(
	core::Size res1,
	core::Size res2
) const
{

	bool res1_in_region1( region1_.find( res1 ) != region1_.end() );
	bool res1_in_region2( region2_.find( res1 ) != region2_.end() );

	bool res2_in_region1( region1_.find( res2 ) != region1_.end() );
	bool res2_in_region2( region2_.find( res2 ) != region2_.end() );


	if ( (res1_in_region1 && res2_in_region2) || (res2_in_region1 && res1_in_region2 ) ) {
		num_cons_++;
	}

}


PackerNeighborGraphFilter::PackerNeighborGraphFilter()
: task_( /* NULL */ ), sfxn_( nullptr ), task_invalidated_( false )
{
	required_connections_per_res_.clear();
	required_connections_between_regions_.clear();
}

PackerNeighborGraphFilter::~PackerNeighborGraphFilter()= default;

void
PackerNeighborGraphFilter::set_required_connections_for_residue(
	core::Size residue,
	core::Size required_connections
){

	std::pair< core::Size, core::Size > con_pair( std::make_pair( residue, required_connections ) );

	required_connections_per_res_.insert( con_pair );
}


void
PackerNeighborGraphFilter::add_required_connection_for_residue(
	core::Size residue
){

	auto map_it = required_connections_per_res_.find( residue );

	if ( map_it == required_connections_per_res_.end() ) {
		required_connections_per_res_.insert( std::make_pair( residue, 1 ) );
	} else map_it->second++;
}


void
PackerNeighborGraphFilter::add_required_connections_between_regions(
	std::set< core::Size > const & region1,
	std::set< core::Size > const & region2,
	core::Size required_connections
){
	required_connections_between_regions_.push_back( RegionalConnections( region1, region2, required_connections ) );
}

bool
PackerNeighborGraphFilter::apply( core::pose::Pose const & pose ) const {

	using namespace utility::graph;

	//if( task_invalidated_ ) utility_exit_with_message("Calling PackerNeighborGraphFilter apply function even though the task has been invalidated");

	std::map< core::Size, core::Size > connections_for_important_res;
	auto regions_end = required_connections_between_regions_.end();

	//reset all the regions
	for ( auto reg_it = required_connections_between_regions_.begin();
			reg_it != regions_end; ++reg_it ) {
		reg_it->reset_num_connections();
	}

	//first, we create the packer neighbor graph
	GraphCOP png = core::pack::create_packer_graph( pose, *sfxn_, task_ );

	//well, that was easy!!
	//now we have to go through the packer neighbor graph and see what kind of connections it has
	for ( Node::EdgeListConstIter edge_it( png->const_edge_list_begin() ), edge_end( png->const_edge_list_end() );
			edge_it != edge_end; ++edge_it ) {


		core::Size res1 = (*edge_it)->get_first_node_ind();
		core::Size res2 = (*edge_it)->get_second_node_ind();

		//check if either of these residues have required connections
		if ( required_connections_per_res_.find( res1 ) != required_connections_per_res_.end() ) {

			auto map_it = connections_for_important_res.find( res1 );
			if ( map_it  == connections_for_important_res.end() ) {

				connections_for_important_res.insert( std::pair< core::Size, core::Size >( res1, 1 ) );
			} else map_it->second++;
		}

		if ( required_connections_per_res_.find( res2 ) != required_connections_per_res_.end() ) {

			auto map_it = connections_for_important_res.find( res2 );
			if ( map_it  == connections_for_important_res.end() ) {

				connections_for_important_res.insert( std::pair< core::Size, core::Size >( res2, 1 ) );
			} else map_it->second++;
		}

		//now check whether these two residues belong to regions that we care about
		for ( auto reg_it = required_connections_between_regions_.begin();
				reg_it != regions_end; ++reg_it ) {
			reg_it->check_if_connected_residues_belong_to_regions( res1, res2 );
		}

	} //iteration over graph edges


	//aight. now let's check whether the connections in the png satisfy the requirements
	//first, let's check if all required residues have a sufficient number of connections
	for ( auto const & required_connections_per_re : required_connections_per_res_ ) {

		std::map< core::Size, core::Size >::const_iterator curmap_it = connections_for_important_res.find( required_connections_per_re.first );

		if ( curmap_it == connections_for_important_res.end() ) return false; //residue was required to have connections but had 0

		else if ( curmap_it->second < required_connections_per_re.second ) return false;  //residue had too little connections

	}

	//second, check whether all the regions have the right number of connections
	for ( auto reg_it = required_connections_between_regions_.begin();
			reg_it != regions_end; ++reg_it ) {
		if ( ! reg_it->enough_connections() ) return false;
	}

	//yay! if we've made it till here, that means all the requirements are met
	return true;

} // apply_filter


} // filters
} // protocols
