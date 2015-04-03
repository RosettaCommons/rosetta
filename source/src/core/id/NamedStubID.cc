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
#include <core/id/NamedStubID.hh>

#include <utility/exit.hh>


// C++ headers
#include <ostream>
#include <string>

#include <utility/vector1.hh>
#include <sstream>


namespace core {
namespace id {

// convienience c'stor if the residue is the same for all atoms
NamedStubID::NamedStubID( std::string const& a1, std::string const& a2, std::string const& a3, core::Size rsd ) :
	center_( a1, rsd ),
	atom1( a1, rsd ),
	atom2( a2, rsd ),
	atom3( a3, rsd )
{}

NamedStubID::NamedStubID( std::string const& c, std::string const& a1, std::string const& a2, std::string const& a3, core::Size rsd ) :
	center_( c, rsd ),
	atom1( a1, rsd ),
	atom2( a2, rsd ),
	atom3( a3, rsd )
{}

// convienience c'stor
NamedStubID::NamedStubID( std::string const& a1, Size rsd1, std::string const& a2, Size rsd2, std::string a3, Size rsd3 ) :
	center_( a1, rsd1 ),
	atom1( a1, rsd1 ),
	atom2( a2, rsd2 ),
	atom3( a3, rsd3 )
{}

NamedStubID::NamedStubID( AtomList const& atoms, core::Size rsd ) {
debug_assert( atoms.size() == 3 || atoms.size() == 4 );
	Size ind;
	if ( atoms.size() == 4 ) {
		center_.atom() = atoms[ 1 ];
		center_.rsd() = rsd;
		ind = 2;
	} else ind = 1;
	atom1.atom() = atoms[ ind++ ];
	atom2.atom() = atoms[ ind++ ];
	atom3.atom() = atoms[ ind++ ];
	atom1.rsd() = rsd;
	atom2.rsd() = rsd;
	atom3.rsd() = rsd;
}

NamedAtomID const &
NamedStubID::atom( Size const index ) const
{
  switch ( index ) {
  case 1:
    return atom1;
  case 2:
    return atom2;
  case 3:
    return atom3;
  default:
    utility_exit_with_message("StubID's have exactly three atoms, 1-3");
  }
  return atom1; // won't get here
}

/// @brief stream << NamedStubID
std::ostream &
operator <<( std::ostream& os, NamedStubID const & a )
{
  os << "STUB: " << a.center_ << ' ' << a.atom1 << ' ' << a.atom2 << ' ' << a.atom3 << ' ';
  return os;
}

/////////////////////////////////////////////////////////////////////////////
/// @brief input operator >> NamedAtomID
std::istream &
operator >>( std::istream & is, NamedStubID& s )
{
	std::string tag;
	is >> tag;
	if ( tag != "STUB:" ) {
		is.setstate( std::ios_base::failbit );
		return is;
	}
	is >> s.center_ >> s.atom1 >> s.atom2 >> s.atom3;
	return is;
}


} // namespace id
} // namespace core
