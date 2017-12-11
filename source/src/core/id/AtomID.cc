// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#endif // SERIALIZATION


namespace core {
namespace id {

static basic::Tracer tr( "core.id.AtomID" );


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


/// @brief stream << AtomID
std::ostream &
operator <<( std::ostream& os, StubID const & stub )
{
	os << " Stub: a1 " << stub.atom1 << " a2 " << stub.atom2 << " a3 " << stub.atom3 << " cen " << stub.center_;
	return os;
}

bool StubID::operator == ( StubID const & rhs ) const
{
	return atom1 == rhs.atom1 && atom2 == rhs.atom2 && atom3 == rhs.atom3 && center_ == rhs.center_;
}

bool StubID::operator != ( StubID const & rhs ) const
{
	return ! ( *this == rhs );
}

} // namespace id
} // namespace core

///@brief Set the value of atom and residue.
void
core::id::AtomID::set( core::Size atomno_in, core::Size rsd_in){
	atomno_ = atomno_in;
	rsd_ = rsd_in;
}


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::id::AtomID::save( Archive & arc ) const {
	arc( CEREAL_NVP( atomno_ ) ); // Size
	arc( CEREAL_NVP( rsd_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::id::AtomID::load( Archive & arc ) {
	arc( atomno_ ); // Size
	arc( rsd_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::id::AtomID );

#endif // SERIALIZATION

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::id::BondID::save( Archive & arc ) const {
	arc( CEREAL_NVP( atom1 ) ); // class core::id::AtomID
	arc( CEREAL_NVP( atom2 ) ); // class core::id::AtomID
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::id::BondID::load( Archive & arc ) {
	arc( atom1 ); // class core::id::AtomID
	arc( atom2 ); // class core::id::AtomID
}

SAVE_AND_LOAD_SERIALIZABLE( core::id::BondID );

#endif // SERIALIZATION

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::id::StubID::save( Archive & arc ) const {
	arc( CEREAL_NVP( atom1 ) ); // class core::id::AtomID
	arc( CEREAL_NVP( atom2 ) ); // class core::id::AtomID
	arc( CEREAL_NVP( atom3 ) ); // class core::id::AtomID
	arc( CEREAL_NVP( center_ ) ); // class core::id::AtomID
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::id::StubID::load( Archive & arc ) {
	arc( atom1 ); // class core::id::AtomID
	arc( atom2 ); // class core::id::AtomID
	arc( atom3 ); // class core::id::AtomID
	arc( center_ ); // class core::id::AtomID
}

SAVE_AND_LOAD_SERIALIZABLE( core::id::StubID );

#endif // SERIALIZATION
