// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/simple_filters/DomainInterfaceFilter.cc
/// @brief Implementation of DomainInterfaceFilter which checks whether a defined region of a protein
/// is part of an interface with a pre-defined other domain of the protein
///
/// @author Robert Lindner <rlindner@mpimf-heidelberg.mpg.de>

// Unit Headers
#include <protocols/simple_filters/DomainInterfaceFilter.hh>

// Package Headers
#include <basic/datacache/DataMap.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/pose/selection.hh>
#include <core/pack/task/operation/util/interface_vector_calculate.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/tag/Tag.hh>

//// C++ headers
#include <string>

namespace protocols {
namespace simple_filters {

DomainInterfaceFilter::DomainInterfaceFilter():
	filters::Filter( "DomainInterfaceFilter" ),
	cb_dist_cut_( 11.0 ),
	nearby_atom_cut_( 5.5 ),
	vector_angle_cut_( 75.0 ),
	vector_dist_cut_( 9.0 )
{ }

DomainInterfaceFilter::DomainInterfaceFilter( DomainInterfaceFilter const & src ):
	filters::Filter( src.type_ ),
	cb_dist_cut_( src.cb_dist_cut_ ),
	nearby_atom_cut_( src.nearby_atom_cut_ ),
	vector_angle_cut_( src.vector_angle_cut_ ),
	vector_dist_cut_( src.vector_dist_cut_ ),
	query_region_( src.query_region_ ),
	target_region_( src.target_region_ ),
	query_region_str_( src.query_region_str_ ),
	target_region_str_( src.target_region_str_ )
{}


DomainInterfaceFilter::~DomainInterfaceFilter() {
}

  
void DomainInterfaceFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &
) {
	cb_dist_cut( tag->getOption< core::Real >( "cb_dist_cut", 11.0 )); // should be bigger than vector_dist_cut
	nearby_atom_cut( tag->getOption< core::Real >( "nearby_atom_cut", 5.5 ) );
	vector_angle_cut( tag->getOption< core::Real >( "vector_angle_cut", 75.0 ) );
	vector_dist_cut( tag->getOption< core::Real >( "vector_dist_cut", 9.0 ) );
	
	if( tag->hasOption( "target" ) ) {
		target_region_string( tag->getOption< std::string >( "target" ) );
	}
	if( tag->hasOption( "query" ) ) {
		query_region_string( tag->getOption< std::string >( "query" ) );
	}

}
	
filters::FilterOP DomainInterfaceFilter::clone() const {
	return new DomainInterfaceFilter( *this );
}

filters::FilterOP DomainInterfaceFilter::fresh_instance() const {
	return new DomainInterfaceFilter();
}


/// @brief Returns true if the given pose passes the filter, false otherwise.
/// Pose passes the filter if ALL query residues are part of an interface with
/// ANY of the target residues.
bool DomainInterfaceFilter::apply( core::pose::Pose const & pose ) const {
	std::set< core::Size > query_local, target_local;

	// set up query set	
	if( !query_region_.empty() ) 
		query_local.insert( query_region_.begin(), query_region_.end() );
	else if( !query_region_str_.empty() ) 
		 query_local = core::pose::get_resnum_list( query_region_str_, pose );

	// set up target set	
	if( !target_region_.empty() ) 
		target_local.insert( target_region_.begin(), target_region_.end() );
	else if( !target_region_str_.empty() ) 
		 target_local = core::pose::get_resnum_list( target_region_str_, pose );

	utility::vector1< bool > rs_subset( false, pose.total_residue() );
	rs_subset = core::pack::task::operation::util::calc_interacting_vector( pose, target_local, query_local, cb_dist_cut_, nearby_atom_cut_, vector_angle_cut_, vector_dist_cut_ );
	
	for( std::set< core::Size >::const_iterator it = query_local.begin();
		it != query_local.end();
		++it)
	{
		if( *it <= rs_subset.size() && !rs_subset[ *it ] ) return false;
	}
	return true;	
}

// Interface Selector Parameters
core::Real DomainInterfaceFilter::cb_dist_cut() const { return cb_dist_cut_; }
core::Real DomainInterfaceFilter::nearby_atom_cut() const { return nearby_atom_cut_; }
core::Real DomainInterfaceFilter::vector_angle_cut() const { return vector_angle_cut_; }
core::Real DomainInterfaceFilter::vector_dist_cut() const { return vector_dist_cut_; }
void DomainInterfaceFilter::cb_dist_cut( core::Real setting ) { cb_dist_cut_ = setting; }
void DomainInterfaceFilter::nearby_atom_cut( core::Real setting ) { nearby_atom_cut_ = setting; }
void DomainInterfaceFilter::vector_angle_cut( core::Real setting ) { vector_angle_cut_ = setting; }
void DomainInterfaceFilter::vector_dist_cut( core::Real setting ) { vector_dist_cut_ = setting; }

// Setters and getters for the query and target strings
std::string const & DomainInterfaceFilter::query_region_string() { return query_region_str_; }
std::string const & DomainInterfaceFilter::target_region_string() { return target_region_str_; }
void DomainInterfaceFilter::query_region_string( std::string const & query_region_str  ) { query_region_str_ = query_region_str; }
void DomainInterfaceFilter::target_region_string( std::string const & target_region_str  ) { target_region_str_ = target_region_str; }

std::set< core::Size > const & DomainInterfaceFilter::query_region() { return query_region_; }
std::set< core::Size > const & DomainInterfaceFilter::target_region() { return target_region_; }
void DomainInterfaceFilter::query_region( std::set< core::Size > & region ) { query_region_ = region; }
void DomainInterfaceFilter::target_region( std::set< core::Size > & region ) { target_region_ = region; }

} // simple_filters
} // protocols

