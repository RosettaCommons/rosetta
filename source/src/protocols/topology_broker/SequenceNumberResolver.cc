// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SequenceNumberResolver
/// @brief
/// @author Justin Porter, Tatjana Braun

// Unit Headers
#include <protocols/topology_broker/SequenceNumberResolver.hh>

// Package Headers

// Project Headers

// ObjexxFCL Headers

// Utility headers
#include <utility/excn/Exceptions.hh>

// utility
#include <basic/Tracer.hh>

// C/C++ headers
#include <sstream>
#include <map>
#include <utility>


#ifdef WIN32
#include <iterator>
#endif

static thread_local basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

SequenceNumberResolver::SequenceNumberResolver( const SequenceNumberResolver&  src ):
	ReferenceCount()
{
	offset_map_ = src.offset_map_;
	offset_map_reversed_ = src.offset_map_reversed_;
}

/// @brief Returns offset of a given label.
core::Size
SequenceNumberResolver::offset( std::string const& label ) const {
	OffsetMap::const_iterator p = offset_map_.find(label);

	if( p != offset_map_.end() ) { return p->second; }
	//else if( label == "" ) {
	//  //Zero offset sequences labels should be "BASE" - an empty label probably means a logic bug somewhere
	//	tr.Warning << "Warning: Attempting to resolve a sequence number with an empty claim label - assuming zero offset." << std::endl;
	//	return 0;
	} else {
		throw utility::excn::EXCN_BadInput( "SequenceNumberResolver asked to resolve SequenceClaim label '"
						+ label + "', which does not match any SequenceClaim labels." );
	}
}

// offset = number of residues in front of sequence claim with specified label
// offset + 1 = global position of first element of sequence claim with specified label
void SequenceNumberResolver::register_label_offset( std::string const& label, core::Size offset ){

	std::pair<std::map< std::string, core::Size >::iterator,bool> ret;
	ret = offset_map_.insert( std::pair < std::string, core::Size >( label,offset ) );
	if ( ret.second == false ){
		throw utility::excn::EXCN_BadInput( "Multiple sequence claims with label " + label );
	}

	//Store entry as well in reversed map (key=offset, value=label)
	std::pair<std::map< core::Size, std::string  >::iterator,bool> ret_rev;
	ret_rev = offset_map_reversed_.insert( std::pair < core::Size, std::string >( offset, label ) );
	if ( ret_rev.second == false ){
		std::ostringstream msg;
		msg << "Multiple sequence claims with offset " << offset;
		throw utility::excn::EXCN_BadInput( msg.str() );
	}
}

core::Size SequenceNumberResolver::find_global_pose_number(std::string const&label, core::Size resid) const {
	return offset(label) + resid ;
}

core::Size SequenceNumberResolver::find_global_pose_number(std::pair< std::string, core::Size> pos_pair ) const {
	return find_global_pose_number( pos_pair.first, pos_pair.second );
}


core::Size SequenceNumberResolver::find_global_pose_number(std::string const&label) const {
	return find_global_pose_number( label, 1 );
}

std::string SequenceNumberResolver::find_label( core::Size pose_number ) const {
	std::map<core::Size, std::string>::const_iterator itlow = search_reversed_map( pose_number );
	return itlow->second;
}

core::Size SequenceNumberResolver::find_local_pose_number( core::Size pose_number ) const{
	std::map<core::Size, std::string>::const_iterator itlow = search_reversed_map( pose_number );
	return ( pose_number - itlow->first );
}

std::pair<core::Size, std::string> SequenceNumberResolver::terminal_pair() const {
	return std::make_pair( offset_map_reversed_.rbegin()->first, offset_map_reversed_.rbegin()->second );
}


std::map<core::Size, std::string>::const_iterator SequenceNumberResolver::search_reversed_map( core::Size pose_number ) const{

	//Returns first entry whose key is greater or equal compared to the given pose_number
	std::map<core::Size, std::string>::const_iterator itlow = offset_map_reversed_.lower_bound( pose_number ) ;

	//Decrement iterator to access element that is smaller than given pose number
	if ( itlow != offset_map_reversed_.begin() ) itlow--;
	else {
		std::ostringstream msg;
		msg << "Iterator of SequenceNumberResolver is out of range for pose number " << pose_number;
		throw utility::excn::EXCN_RangeError( msg.str() );
	}
	return itlow;
}

// void SequenceNumberResolver::clear(){
// 	offset_map_.clear();
// 	offset_map_reversed_.clear();
// }

}
}
