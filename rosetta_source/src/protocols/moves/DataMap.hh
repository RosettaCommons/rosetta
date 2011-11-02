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

#ifndef INCLUDED_protocols_moves_DataMap_hh
#define INCLUDED_protocols_moves_DataMap_hh

// Unit Headers
// AUTO-REMOVED #include <protocols/moves/Mover.fwd.hh>

// Package headers
// AUTO-REMOVED #include <protocols/moves/MoverStatistics.hh>
// AUTO-REMOVED #include <protocols/moves/MoverStatus.hh>

// Project headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <utility/tag/Tag.fwd.hh>
// AUTO-REMOVED #include <utility/exit.hh>
#include <sstream>
// ObjexxFCL Headers

// C++ Headers
// AUTO-REMOVED #include <string>
#include <map>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/exit.hh>
#include <string>

namespace protocols {
namespace moves {

/// @brief general-purpose store for any reference-count derived object

class DataMap : public utility::pointer::ReferenceCount {
public:
	typedef std::map< std::string, std::map< std::string, utility::pointer::ReferenceCountOP > >::iterator iterator;
	typedef std::map< std::string, std::map< std::string, utility::pointer::ReferenceCountOP > >::const_iterator const_iterator;

public:
	DataMap();
	~DataMap();
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
template< class Ty >
Ty
DataMap::get( std::string const type, std::string const name ) const {
	using namespace utility::pointer;
	Ty ret( 0 );
	
	if( !has( type, name ) ){
		std::stringstream error_message;
		error_message << "ERROR: Could not find "<<type<<" and name "<<name<<" in Datamap\n";
		
		utility_exit_with_message( error_message.str() );
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
		
		utility_exit_with_message( error_message.str() );
	}
	return( ret );
}

} // moves
} // protocols

#endif
