// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/ContactTypes.cc
/// @brief  fragment picker contact types
/// @author David E. Kim (dekim@u.washington.edu)

// package headers
#include <protocols/frag_picker/ContactTypes.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>

namespace protocols {
namespace frag_picker {


/// BEGIN: local functions

/// @brief setup the map that converts string name to enum
std::map< std::string, ContactType > setup_name2type() {
	std::map< std::string, ContactType > n2t;
	n2t[ "ca" ] = CA;
	n2t[ "cb" ] = CB;
	n2t[ "cen" ] = CEN;
	n2t[ "unk" ] = UNK;
	return n2t;
}

/// @brief map that converts string name to enum
inline
std::map< std::string, ContactType > & name2type() {
	// static initialization only happens once
	static std::map< std::string, ContactType > * name2contacttype_ = new std::map< std::string, ContactType >( setup_name2type() );
	return *name2contacttype_;
}

/// @brief setup the vector that maps enum to string name
utility::vector1< std::string > setup_type2name() {
	utility::vector1< std::string > t2n( num_contact_types );
	for ( std::map< std::string, ContactType >::const_iterator iter = name2type().begin(),
			iter_end = name2type().end(); iter != iter_end; ++iter ) {
		t2n[ iter->second ] = iter->first;
	}
	return t2n;
}

/// @brief vector that maps enum to string name
inline
utility::vector1< std::string > & type2name() {
	// static initialization only happens once
	static utility::vector1< std::string > * contacttype2name_ = new utility::vector1< std::string >( setup_type2name() );
	return *contacttype2name_;
}

/// END: local functions

//////////////////////////////////////////////////////////
/// @brief give a string name and return its enum type
//////////////////////////////////////////////////////////
ContactType
contact_type( std::string const & name ) {
	std::map< std::string, ContactType >::const_iterator iter = name2type().find( name );
	if ( iter == name2type().end() ) {
		utility_exit_with_message( "unrecognized contact type " + name );
	}
	return iter->second;
}

///////////////////////////////////////////////////////
/// @brief give an enum type and return the string name
///////////////////////////////////////////////////////
std::string
contact_name( ContactType type ) {
	if ( type > num_contact_types ) return "ContactTypeOutofRange";
	return type2name()[ type ];
}


} // namespace frag_picker
} // namespace protocols

