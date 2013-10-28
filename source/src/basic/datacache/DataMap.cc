// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Sarel Fleishman

// Unit Headers
#include <basic/datacache/DataMap.hh>

// Package headers

// Project headers
// ObjexxFCL Headers

// C++ Headers
// AUTO-REMOVED #include <string>
#include <map>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <basic/Tracer.hh>

namespace basic {
namespace datacache {

static basic::Tracer TR( "basic.datacache.DataMap" );

DataMap::DataMap() {}
DataMap::~DataMap() {}

bool
DataMap::add( std::string const type, std::string const name, utility::pointer::ReferenceCountOP const op ){
	if( has( type, name ) ){
		TR<<"A datum of type "<<type<<" and name "<<name<<" has been added before. I'm not adding again. This is probably a BIG error but I'm letting it pass!"<<std::endl;
		return false;
	}

	data_map_[ type ].insert(std::make_pair( name, op ) );
	return true;
}

bool
DataMap::has( std::string const type, std::string const name/*=""*/ ) const {
	std::map< std::string, std::map< std::string, utility::pointer::ReferenceCountOP > >::const_iterator it;

	it = data_map_.find( type );
	if( it == data_map_.end() ) return false;
	std::map< std::string, utility::pointer::ReferenceCountOP >::const_iterator it2;
	it2 = it->second.find( name );
	if( it2 == it->second.end() ) return false;

	return true;
}

std::map< std::string, utility::pointer::ReferenceCountOP > &
DataMap::operator []( std::string const & type ) {
	if( !has( type ) ) {
// "dummy_entry" serves as a placeholder while the datamap does not contain actual maps of this type.
// it is removed if the map is accessed.
		add( type, "dummy_entry", 0 );
	}

	std::map< std::string, utility::pointer::ReferenceCountOP > & m( data_map_.find( type )->second );
	if ( m.size() > 1 ) {
		std::map< std::string, utility::pointer::ReferenceCountOP >::iterator it;
		it=m.find( "dummy_entry" );
		m.erase( it );
	}
	return m;
}

platform::Size
DataMap::size() const { return data_map_.size(); }

DataMap::iterator
DataMap::begin() { return data_map_.begin(); }

DataMap::iterator
DataMap::end() { return data_map_.end(); }

DataMap::const_iterator
DataMap::begin() const { return data_map_.begin(); }

DataMap::const_iterator
DataMap::end() const { return data_map_.end(); }

} // datacache
} // basic
