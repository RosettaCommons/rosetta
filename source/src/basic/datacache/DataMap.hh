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

#ifndef INCLUDED_basic_datacache_DataMap_hh
#define INCLUDED_basic_datacache_DataMap_hh

// Project headers
#include <core/types.hh>
#include <sstream>
// ObjexxFCL Headers

// C++ Headers
// AUTO-REMOVED #include <string>
#include <map>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

//#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <string>
#include <basic/Tracer.hh>

namespace basic {
namespace datacache {

static basic::Tracer TR_hh( "basic.datacache.DataMap_hh" );

/// @brief general-purpose store for any reference-count derived object

class DataMap : public utility::pointer::ReferenceCount {
public:
	typedef std::map< std::string, std::map< std::string, utility::pointer::ReferenceCountOP > >::iterator iterator;
	typedef std::map< std::string, std::map< std::string, utility::pointer::ReferenceCountOP > >::const_iterator const_iterator;

public:
	DataMap();
	virtual ~DataMap();
	iterator begin();
	iterator end();
	const_iterator begin() const;
	const_iterator end() const;
	// @brief add an object to the map, returning false if an object of that
	// name already exists
	bool add(
		std::string const type,
		std::string const name,
		utility::pointer::ReferenceCountOP const op
	);
	bool has( std::string const type, std::string const name="" ) const;
	template< class Ty > Ty get( std::string const type, std::string const name ) const;
	std::map< std::string, utility::pointer::ReferenceCountOP > & operator [](
		std::string const & type
	);
	/// @brief returns the size of the map (how many different types are in data_map_
	core::Size size() const;

private:
	std::map< std::string, std::map< std::string, utility::pointer::ReferenceCountOP > > data_map_;
};

/// @details a template utility function to grab any type of object from the
/// Data_map. Downcasts the ReferenceCount object in map to the template data
/// type using dynamic_cast to ensure type-correctness
/// @throws Throws a utility::excn::EXCN_Msg_Exception in the event that
/// the requested object cannot be found in the DataMap.
template< class Ty >
Ty
DataMap::get( std::string const type, std::string const name ) const {
	using namespace utility::pointer;
	Ty ret( 0 );

	if( !has( type, name ) ){
		std::stringstream error_message;
		error_message << "ERROR: Could not find "<<type<<" and name "<<name<<" in Datamap\n";
		throw utility::excn::EXCN_Msg_Exception( error_message.str() );
	}

	std::map< std::string, utility::pointer::ReferenceCountOP > const dm( data_map_.find( type )->second );
	for( std::map< std::string, ReferenceCountOP >::const_iterator it=dm.begin(); it!=dm.end(); ++it ) {
		if( it->first == name ) {
			ret = dynamic_cast< Ty >( it->second.get() );
			break;
		}
	}
	if( ret==0 ) {
		std::stringstream error_message;
		error_message << "ERROR: Dynamic_cast failed for type "<<type<<" and name "<<name<<'\n';
		throw utility::excn::EXCN_Msg_Exception( error_message.str() );
	}
	return( ret );
}


/// @brief templated function for adding or getting an item from the datamap. Automatically checks whether an item
/// of the requested type and name exists on the datamap. If so, returns the OP for that item, if not, instantiates
/// that item on the datamap and returns the OP for it.
template < class Ty >
Ty *
get_set_from_datamap( std::string const type, std::string const name, basic::datacache::DataMap & data ){
	Ty *obj;
	if( data.has( type, name ) ){
		obj = data.get< Ty * >( type, name );
		TR_hh<<"Getting object-type, name "<<type<<' '<<name<<" from datamap"<<std::endl;
	}
	else{
		obj = new Ty;
		data.add( type, name, obj );
		TR_hh<<"Adding object-type, name "<<type<<' '<<name<<" to datamap"<<std::endl;
	}
	return obj;
}

} // datacache
} // basic

#endif
