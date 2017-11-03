// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

static THREAD_LOCAL basic::Tracer tr( "core.id.NamedAtomID", basic::t_info );

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#endif // SERIALIZATION


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


} // namespace id
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::id::NamedAtomID::save( Archive & arc ) const {
	arc( CEREAL_NVP( atom_ ) ); // std::string
	arc( CEREAL_NVP( rsd_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::id::NamedAtomID::load( Archive & arc ) {
	arc( atom_ ); // std::string
	arc( rsd_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::id::NamedAtomID );

#endif // SERIALIZATION
