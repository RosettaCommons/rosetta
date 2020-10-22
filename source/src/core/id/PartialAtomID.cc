// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/PartialAtomID.cc
/// @brief  Partially-resolved atom identifier
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit headers
#include <core/id/PartialAtomID.hh>
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

static basic::Tracer tr( "core.id.PartialAtomID" );

///@brief Set the value of atom and residue.
void
PartialAtomID::set_complete( core::Size atomno_in, core::Size rsd_in){
	atomno_ = atomno_in;
	resconnid_ = 0;
	bonds_from_resconn_ = 0;
	rsd_ = rsd_in;
}

void
PartialAtomID::set_partial( core::Size resconnid_in, core::Size bonds_from_resconn_in, core::Size rsd_in){
	atomno_ = 0;
	resconnid_ = resconnid_in;
	bonds_from_resconn_ = bonds_from_resconn_in;
	rsd_ = rsd_in;
}


/// @brief stream << PartialAtomID
std::ostream &
operator <<( std::ostream& os, PartialAtomID const & a )
{
	os << "atomno= " << a.atomno() << " resconn_id= " << a.resconnid() << " bonds_from_resconn= " << a.bonds_from_resconn() << " rsd= " << a.rsd();
	return os;
}


} // namespace id
} // namespace core



#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::id::PartialAtomID::save( Archive & arc ) const {
	arc( CEREAL_NVP( atomno_ ) ); // Size
	arc( CEREAL_NVP( resconnid_ ) ); // Size
	arc( CEREAL_NVP( bonds_from_resconn_ ) ); // Size
	arc( CEREAL_NVP( rsd_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::id::PartialAtomID::load( Archive & arc ) {
	arc( atomno_ ); // Size
	arc( resconnid_ ); // Size
	arc( bonds_from_resconn_ ); // Size
	arc( rsd_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::id::PartialAtomID );

#endif // SERIALIZATION

