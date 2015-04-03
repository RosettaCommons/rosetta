// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/AtomID.cc
/// @brief  Kinematics Atom identifier
/// @author Phil Bradley


// Unit headers
#include <core/id/AtomID.hh>
//#include <core/chemical/ResidueType.hh>

#include <utility/exit.hh>


#include <basic/Tracer.hh>

// C++ headers
#include <ostream>

#include <utility/vector1.hh>


namespace core {
namespace id {

static thread_local basic::Tracer tr( "core.id.AtomID" );


AtomID const &
BondID::other_atom( AtomID const & id ) const
{
	if ( id == atom1 ) return atom2;
	else if ( id == atom2 ) return atom1;
	else utility_exit_with_message( "BondID::other_atom: unknown atom" );
	return atom1; // wont get here
}

/// @brief stream << AtomID
std::ostream &
operator <<( std::ostream& os, AtomID const & a )
{
	os << " atomno= " << a.atomno() << " rsd= " << a.rsd() << ' ';
	return os;
}

/// @brief Globals
AtomID const BOGUS_ATOM_ID( 0,0 );
AtomID const CHAINBREAK_BOGUS_ATOM_ID( 0,0 );
StubID const BOGUS_STUB_ID( BOGUS_ATOM_ID, BOGUS_ATOM_ID, BOGUS_ATOM_ID );

/// @brief stream << AtomID
std::ostream &
operator <<( std::ostream& os, StubID const & stub )
{
	os << " Stub: a1 " << stub.atom1 << " a2 " << stub.atom2 << " a3 " << stub.atom3 << " cen " << stub.center_;
	return os;
}

} // namespace id
} // namespace core
