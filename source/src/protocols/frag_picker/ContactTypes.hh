// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/ContactTypes.hh
/// @brief  fragment picker type declarations
/// @author David E. Kim (dekim@u.washington.edu)

#ifndef INCLUDED_protocols_frag_picker_ContactTypes_hh
#define INCLUDED_protocols_frag_picker_ContactTypes_hh

// C++ headers
#include <string>

namespace protocols {
namespace frag_picker {

enum ContactType {
	CA=1,
	CB,
	CEN,
	UNK,
	num_contact_types = UNK //keep this guy last
};

//////////////////////////////////////////////////////////
/// @brief give a string name and return its enum type
//////////////////////////////////////////////////////////
ContactType
contact_type( std::string const & name );

///////////////////////////////////////////////////////
/// @brief give an enum type and return the string name
///////////////////////////////////////////////////////
std::string
contact_name( ContactType type );


} // namespace frag_picker
} // namespace protocols

#endif // INCLUDED_protocols_frag_picker_contact_types_HH
