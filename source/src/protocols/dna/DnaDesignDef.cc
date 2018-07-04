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
	// We have to catch this here rather than doing a aa_from_name call later in the taskop.
	// Why? Because we sometimes feed the taskop true name3s, and aa_from_name can't handle
	// that for obvious reasons. This is not fast, but it's outside of a tight loop and
	// is in a function already doing string operations. Perhaps faster than doing the
	// equivalent in the task.
	if ( name3 == "ADE" ) { name3 = " DA"; }
	else if ( name3 == "CYT" ) { name3 = " DC"; }
	else if ( name3 == "GUA" ) { name3 = " DG"; }
	else if ( name3 == "THY" ) { name3 = " DT"; }
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
