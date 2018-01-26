// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Sarel Fleishman

#ifndef INCLUDED_basic_datacache_DataMap_hh
#define INCLUDED_basic_datacache_DataMap_hh

// Project headers
#include <basic/Tracer.hh>                            // for Tracer
#include <map>                                        // for map, __map_cons...
#include <memory>                                     // for shared_ptr, dyn...
#include <platform/types.hh>                          // for Size
#include <sstream>                                    // for string, operator<<
#include <string>                                     // for char_traits
#include <utility/excn/Exceptions.hh>                 // for EXCN_Msg_Exception
#include <utility>                                    // for pair
#include <utility/pointer/owning_ptr.hh>              // for dynamic_pointer_cast
#include <utility/pointer/ReferenceCount.hh>          // for ReferenceCount

namespace basic {
namespace datacache {

static basic::Tracer TR_hh( "basic.datacache.DataMap_hh" );

/// @brief general-purpose store for any reference-count derived object

class DataMap : public utility::pointer::ReferenceCount {
public:
	typedef std::map< std::string, std::map< std::string, utility::pointer::ReferenceCountOP > >::iterator iterator;
	typedef std::map< std::string, std::map< std::string, utility::pointer::ReferenceCountOP > >::const_iterator const_iterator;
	typedef std::map< std::string, utility::pointer::ReferenceCountCOP >::const_iterator resource_const_iterator;

public:
	DataMap();
	~DataMap() override;

	iterator begin();
	iterator end();

	const_iterator begin() const;
	const_iterator end() const;

	resource_const_iterator resources_begin() const;
	resource_const_iterator resources_end() const;

	// @brief add an object to the map, returning false if an object of that
	// name already exists
	virtual bool add(
		std::string const & type,
		std::string const & name,
		utility::pointer::ReferenceCountOP op
	);

	// @brief add a resource to the data map, returning false if a resource of that
	// name already exists
	bool add_resource(
		std::string const & resource_name,
		utility::pointer::ReferenceCountCOP op
	);


	/// @brief Does the data map contain the given type?
	bool has_type( std::string const & type ) const;

	/// @brief Does the data map contain a resource with the given name?
	bool has_resource( std::string const & resource_name ) const;

	/// @brief Does the data map contain an entry with a specific name in the given type?
	/// @note calling this function without providing a name is just plain wrong and makes no sense
	/// and I would be changing that right now if Kale hadn't already found out that some code
	/// relies on this bad behavior back in pull request #187
	bool has( std::string const & type, std::string const & name="" ) const;
	template< class T > T get( std::string const & type, std::string const & name ) const;
	template< class T > utility::pointer::shared_ptr< T > get_ptr( std::string const & type, std::string const & name ) const;
	template< class T > utility::pointer::shared_ptr< T const > get_resource( std::string const & resource_name ) const;

	std::map< std::string, utility::pointer::ReferenceCountOP > & operator [](
		std::string const & type
	);

	std::map< std::string, utility::pointer::ReferenceCountOP > const &
	category_map(
		std::string const & type
	) const;

	/// @brief returns the size of the map (how many different types are in data_map_
	platform::Size size() const;

private:
	std::map< std::string, std::map< std::string, utility::pointer::ReferenceCountOP > > data_map_;
	std::map< std::string, utility::pointer::ReferenceCountCOP > resource_map_;

};

/// @details a template utility function to grab any type of object from the
/// Data_map. Downcasts the ReferenceCount object in map to the template data
/// type using dynamic_cast to ensure type-correctness
/// @throws Throws a utility::excn::EXCN_Msg_Exception in the event that
/// the requested object cannot be found in the DataMap.
template< class T >
T
DataMap::get( std::string const & type, std::string const & name ) const {
	using namespace utility::pointer;
	T ret( 0 );

	if ( !has( type, name ) ) {
		std::stringstream error_message;
		error_message << "ERROR: Could not find "<<type<<" and name "<<name<<" in Datamap\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
	}

	std::map< std::string, utility::pointer::ReferenceCountOP > const & dm( data_map_.find( type )->second );
	auto iter = dm.find( name );

	if ( iter != dm.end() ) {
		ret = dynamic_cast< T >( iter->second.get() );
		if ( ! ret ) {
			std::stringstream error_message;
			error_message << "ERROR: Dynamic_cast failed for type "<<type<<" and name "<<name<<'\n';
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
		}
	} else {
		std::stringstream error_message;
		error_message << "ERROR: Unable to find data with type \""<<type<<"\" and name \""<<name<< "\" in the data map\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
	}
	return ret;
}

/// @details a template utility function to grab any type of object from the
/// Data_map. Downcasts the owning pointer in map to the template data
/// type using dynamic_pointer_cast to ensure type-correctness
/// @throws Throws a utility::excn::EXCN_Msg_Exception in the event that
/// the requested object cannot be found in the DataMap.
template< class T >
utility::pointer::shared_ptr< T >
DataMap::get_ptr( std::string const & type, std::string const & name ) const {
	using namespace utility::pointer;
	utility::pointer::shared_ptr< T > ret( 0 );

	if ( !has( type, name ) ) {
		std::stringstream error_message;
		error_message << "ERROR: Could not find "<<type<<" and name "<<name<<" in Datamap\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
	}

	std::map< std::string, utility::pointer::ReferenceCountOP > const & dm( data_map_.find( type )->second );
	auto iter = dm.find( name );
	if ( iter != dm.end() ) {
		ret = utility::pointer::dynamic_pointer_cast< T >( iter->second );
		if ( ! ret ) {
			std::stringstream error_message;
			error_message << "ERROR: Dynamic_cast failed for type "<<type<<" and name "<<name<<'\n';
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
		}
	} else {
		std::stringstream error_message;
		error_message << "ERROR: Unable to find data with type \""<<type<<"\" and name \""<<name<< "\" in the data map\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
	}
	return ret;
}

template< class T >
utility::pointer::shared_ptr< T const >
DataMap::get_resource( std::string const & resource_name ) const {
	using namespace utility::pointer;
	utility::pointer::shared_ptr< T const > ret( 0 );

	auto iter = resource_map_.find( resource_name );
	if ( iter == resource_map_.end() ) {
		std::stringstream error_message;
		error_message << "ERROR: Could not find a resource with the name \""<< resource_name <<"\" in the Datamap\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
	}

	ret = utility::pointer::dynamic_pointer_cast< T const >( iter->second );

	if ( ret==0 ) {
		std::stringstream error_message;
		error_message << "ERROR: dynamic_cast failed for resource with name "<< resource_name <<"; the actual resource has a different type.\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
	}
	return ret;
}


/// @brief templated function for adding or getting an item from the datamap. Automatically checks whether an item
/// of the requested type and name exists on the datamap. If so, returns the OP for that item, if not, instantiates
/// that item on the datamap and returns the OP for it.
template < class T >
utility::pointer::shared_ptr< T >
get_set_from_datamap( std::string const & type, std::string const & name, basic::datacache::DataMap & data ){
	utility::pointer::shared_ptr< T > obj;
	if ( data.has( type, name ) ) {
		obj = data.get_ptr< T >( type, name );
		TR_hh<<"Getting object-type, name "<<type<<' '<<name<<" from datamap"<<std::endl;
	} else {
		obj = utility::pointer::shared_ptr< T >( new T );
		data.add( type, name, obj );
		TR_hh<<"Adding object-type, name "<<type<<' '<<name<<" to datamap"<<std::endl;
	}
	return obj;
}

} // datacache
} // basic

#endif
