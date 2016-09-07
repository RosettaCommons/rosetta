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
/// @author ashworth

#include <protocols/dna/DnaDesignDef.hh>

#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
using utility::string_split;

#include <sstream>
#include <string>

namespace protocols {
namespace dna {

using namespace core;

DnaDesignDef::~DnaDesignDef()= default;

DnaDesignDef::DnaDesignDef( std::string const & strdef )
: utility::pointer::ReferenceCount()
{
	// command-line dna_defs are of the format "C.501.ADE"
	// split on '.'
	utility::vector1< std::string > parts( string_split( strdef, '.' ) );
	utility::vector1< std::string >::const_iterator part( parts.begin() );
	chain = (*(part++))[0];
	std::istringstream inum_stream( *(part++) );
	inum_stream >> pdbpos;
	if ( part == parts.end() ) return;
	name3 = *part;
}

std::ostream & operator << ( std::ostream & os, DnaDesignDef const & def )
{
	os << def.chain << "." << def.pdbpos << "." << def.name3;
	return os;
}

std::ostream & operator << ( std::ostream & os, DnaDesignDefs const & defs )
{
	for ( auto def( defs.begin() ), end( defs.end() );
			def != end; ++def ) {
		if ( def != defs.begin() ) os << " ";
		os << *def;
	}
	return os;
}

std::ostream & operator << ( std::ostream & os, DnaDesignDefOPs const & defs )
{
	for ( auto def( defs.begin() ), end( defs.end() );
			def != end; ++def ) {
		if ( def != defs.begin() ) os << " ";
		os << **def;
	}
	return os;
}

} // namespace dna
} // namespace protocols
