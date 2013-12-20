// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/simple_filters/BuriedRegionsFilter.cc
/// @brief Implementation of BuriedRegionsFilter which checks whether a defined region of a protein
/// is buried in an interface by evaluating the number of non-water neighbors of each residue
///
/// @author Robert Lindner <rlindner@mpimf-heidelberg.mpg.de>

// Unit Headers
#include <protocols/simple_filters/BuriedRegionsFilter.hh>

// Package Headers
#include <basic/datacache/DataMap.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/pose/selection.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/tag/Tag.hh>

//// C++ headers
#include <string>

namespace protocols {
namespace simple_filters {

BuriedRegionsFilter::BuriedRegionsFilter():
	filters::Filter( "BuriedRegionsFilter" ),
	distance_cutoff_( 6.0 ),
	neighbor_cutoff_( 10 )
{ }

BuriedRegionsFilter::BuriedRegionsFilter( core::Real distance_cutoff, core::Size neighbor_cutoff ) :
	filters::Filter( "BuriedRegionsFilter" ),
	distance_cutoff_( distance_cutoff ),
	neighbor_cutoff_( neighbor_cutoff ) { }

BuriedRegionsFilter::BuriedRegionsFilter( BuriedRegionsFilter const & src ):
	filters::Filter( src.type_ ),
	region_( src.region_ ),
	region_str_( src.region_str_ ),
	distance_cutoff_( src.distance_cutoff_ ),
	neighbor_cutoff_( src.neighbor_cutoff_ ) {}

BuriedRegionsFilter::~BuriedRegionsFilter() {
}

//virtual void BuriedRegionsFilter::clear() {
//		  *this = BuriedRegionsFilter();
//}
  
void BuriedRegionsFilter::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &
) {
}
	
filters::FilterOP BuriedRegionsFilter::clone() const {
	return new BuriedRegionsFilter( *this );
}

filters::FilterOP BuriedRegionsFilter::fresh_instance() const {
	return new BuriedRegionsFilter();
}


/// @brief Returns true if the given pose passes the filter, false otherwise.
bool BuriedRegionsFilter::apply( core::pose::Pose const & pose ) const {
	std::set< core::Size > res_local( core::pose::get_resnum_list( region_str_, pose ) );
	res_local.insert( region_.begin(), region_.end() );

	core::conformation::PointGraphOP pg = new core::conformation::PointGraph;
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg );
	core::conformation::find_neighbors( pg, distance_cutoff_ );

	utility::vector1< Size > non_water_neighbor_count( pose.total_residue(), 0 );

	for(std::set< core::Size >::const_iterator it = res_local.begin();
		it != res_local.end();
		++it ) 
		{
			if( *it > pose.total_residue() || pose.residue(*it).aa() == core::chemical::aa_h2o ) continue; 

			for ( core::conformation::PointGraph::VertexClass::UpperEdgeListIter
					iter = pg->get_vertex(*it).upper_edge_list_begin(),
					iter_end = pg->get_vertex(*it).upper_edge_list_end();
					iter != iter_end; ++iter ) {
				if ( pose.residue( iter->upper_vertex() ).aa() != core::chemical::aa_h2o ) {
					++non_water_neighbor_count[ *it ];
					++non_water_neighbor_count[ iter->upper_vertex() ];
				}
			}
			std::cerr << "BRFilter" << non_water_neighbor_count[ *it ] << "@" << pose.residue( *it ).name3() << *it << std::endl;
			if ( non_water_neighbor_count[ *it ] < neighbor_cutoff_ ) {
				return false;
			}

	}
	// none of the residues we're interested in evaluate to false, hence:
	return true; 
}

core::Real BuriedRegionsFilter::distance_cutoff() { return distance_cutoff_; }
void BuriedRegionsFilter::distance_cutoff( core::Real distance_cutoff ) { distance_cutoff_ = distance_cutoff; }
core::Size BuriedRegionsFilter::neighbor_cutoff() { return neighbor_cutoff_; }
void BuriedRegionsFilter::neighbor_cutoff( core::Size neighbor_cutoff ) { neighbor_cutoff_ = neighbor_cutoff; }
std::string const & BuriedRegionsFilter::region_string() { return region_str_; }
void BuriedRegionsFilter::region_string( std::string const & region_str  ) { region_str_ = region_str; }
std::set< core::Size > const & BuriedRegionsFilter::region() { return region_; }
void BuriedRegionsFilter::region( std::set< core::Size > & region ) { region_ = region; }


} // simple_filters
} // protocols

