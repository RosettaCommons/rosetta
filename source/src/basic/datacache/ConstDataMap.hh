// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/basic/datacache/ConstDataMap.hh
/// @brief  Class for holding constant owning pointers to data so that it can be safely
///         shared between threads without the need for deep-copy semantics.
/// @author Sarel Fleishman
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_basic_datacache_ConstDataMap_hh
#define INCLUDED_basic_datacache_ConstDataMap_hh

// Project headers
#include <platform/types.hh>
#include <sstream>
// ObjexxFCL Headers

// C++ Headers
// AUTO-REMOVED #include <string>
#include <map>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

//#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <string>
#include <basic/Tracer.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {

/// @brief general-purpose store for any kind of object with the particular
/// copy semantics of copying by value.  This is effectively a map of
/// string pairs to (constant) pointers.  The first string represents the
/// category of the object, and the second being a name for that particular
/// object.  The guarantee with the ConstDataMap is that if an object is put
/// into the map, it may be read from, but it will not be changed underneath
/// you.  Data stored in the ConstDataMap can safely be shared between threads.
class ConstDataMap : public utility::pointer::ReferenceCount {
public:
	typedef std::map< std::string, utility::pointer::ReferenceCountCOP > NamedConstObjectMap;
	typedef std::map< std::string, NamedConstObjectMap > CategorizedConstObjectMap;
	typedef CategorizedConstObjectMap::iterator iterator;
	typedef CategorizedConstObjectMap::const_iterator const_iterator;

public:

	ConstDataMap();
	ConstDataMap( ConstDataMap const & src );
	~ConstDataMap() override;

	/// @brief Performs a shallow copy of all of the pointers stored in rhs into this.
	ConstDataMap & operator = ( ConstDataMap const & rhs );

	/// @brief Performs pointer comparison to determine if these two maps point at the same data.
	bool operator== ( ConstDataMap const & rhs ) const;


	iterator begin();
	iterator end();

	const_iterator begin() const;
	const_iterator end() const;

	// @brief add an object to the map; this will overwrite any existing object if it is present
	void add(
		std::string const & type,
		std::string const & name,
		utility::pointer::ReferenceCountCOP const op
	);

	/// @brief are there any objects in the outer map with the given category?
	bool has( std::string const & category ) const;

	/// @brief Is there an object with the given category and the given name?
	bool has( std::string const & category, std::string const & name ) const;

	template< class Ty >
	Ty const &
	get( std::string const & type, std::string const & name ) const;

	template< class Ty >
	utility::pointer::shared_ptr< Ty const >
	get_ptr( std::string const & type, std::string const & name ) const;

	NamedConstObjectMap &
	operator [](
		std::string const & type
	);

	/// @brief returns the number of objects contained in the map
	platform::Size size() const;

private:
	CategorizedConstObjectMap data_map_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/// @details a template utility function to grab any type of object from the
/// Data_map. Downcasts the ReferenceCount object in map to the template data
/// type using dynamic_cast to ensure type-correctness
/// @throws Throws a utility::excn::EXCN_Msg_Exception in the event that
/// the requested object cannot be found in the ConstDataMap.
template< class Ty >
Ty const &
ConstDataMap::get( std::string const & type, std::string const & name ) const {
	using namespace utility::pointer;

	if ( !has( type, name ) ) {
		std::stringstream error_message;
		error_message << "ERROR: Could not find "<<type<<" and name "<<name<<" in Datamap\n";
		throw utility::excn::EXCN_Msg_Exception( error_message.str() );
	}

	NamedConstObjectMap const & dm( data_map_.find( type )->second );
	auto it = dm.find( name );
	utility::pointer::shared_ptr< Ty const > data_ptr( dynamic_pointer_cast< Ty const > ( it->second ) );
	if ( ! data_ptr ) {
		std::stringstream error_message;
		error_message << "Dynamic_cast failed for type " << type << " and name "<<name<<'\n';
		throw utility::excn::EXCN_Msg_Exception( error_message.str() );
	}
	return *data_ptr;
}

/// @details a template utility function to grab any type of object from the
/// ConstDataMap. Downcasts the owning pointer in map to the template data
/// type using dynamic_pointer_cast to ensure type-correctness
/// @throws Throws a utility::excn::EXCN_Msg_Exception in the event that
/// the requested object cannot be found in the ConstDataMap.
template< class Ty >
utility::pointer::shared_ptr< Ty const >
ConstDataMap::get_ptr( std::string const & type, std::string const & name ) const {
	using namespace utility::pointer;
	utility::pointer::shared_ptr< Ty > ret( 0 );

	if ( !has( type, name ) ) {
		std::stringstream error_message;
		error_message << "ERROR: Could not find "<<type<<" and name "<<name<<" in Datamap\n";
		throw utility::excn::EXCN_Msg_Exception( error_message.str() );
	}

	NamedConstObjectMap const & dm( data_map_.find( type )->second );
	auto it = dm.find( name );
	utility::pointer::shared_ptr< Ty const > data_ptr( dynamic_pointer_cast< Ty const > ( it->second ) );
	if ( ! data_ptr ) {
		std::stringstream error_message;
		error_message << "Dynamic_cast failed for type " << type << " and name "<<name<<'\n';
		throw utility::excn::EXCN_Msg_Exception( error_message.str() );
	}

	return data_ptr;
}

} // datacache
} // basic

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( basic_datacache_ConstDataMap )
#endif // SERIALIZATION


#endif
