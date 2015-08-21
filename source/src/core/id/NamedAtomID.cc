// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/NamedAtomID.cc
/// @brief  Kinematics Atom identifier
/// @author Phil Bradley


// Unit headers
#include <core/id/NamedAtomID.hh>


#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// C++ headers
#include <ostream>

static thread_local basic::Tracer tr( "core.id.NamedAtomID", basic::t_info );

namespace core {
namespace id {

std::string NamedAtomID::to_string() const {
	return atom()+" "+ObjexxFCL::string_of( rsd() );
}

/// @brief stream << NamedAtomID
std::ostream &
operator <<( std::ostream& os, NamedAtomID const & a )
{
	// os << " atom= " << a.atom() << " rsd= " << a.rsd() << ' ';
	os << ObjexxFCL::format::A(4,a.atom()) << " " << ObjexxFCL::format::RJ(4,a.rsd());
	return os;
}


/////////////////////////////////////////////////////////////////////////////
/// @brief input operator >> NamedAtomID
std::istream &
operator >>( std::istream & is, NamedAtomID& e )
{
	std::string tag;
	is >> tag;
	if ( tag != "atom=" ) {
		is.setstate( std::ios_base::failbit );
		return is;
	}
	is >> e.atom();
	is >> tag;
	if ( tag != "rsd=" ) {
		is.setstate( std::ios_base::failbit );
		return is;
	}
	is >> e.rsd();
	return is;
}

/// @brief Globals
NamedAtomID const BOGUS_NAMED_ATOM_ID( "", 0 );
NamedAtomID const CHAINBREAK_BOGUS_NAMED_ATOM_ID( "", 0 );


} // namespace id
} // namespace core
