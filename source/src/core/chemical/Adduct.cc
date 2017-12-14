// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Adduct.cc
/// @brief definition of residue adduct class
/// @author Jim Havranek


// Unit headers
#include <core/chemical/Adduct.hh>

// Utility headers
#include <utility>
#include <utility/exit.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

Adduct::Adduct():
	adduct_name_(),
	atom_name_(),
	atom_type_name_(),
	//atom_type_(),
	mm_atom_type_name_(),
	//mm_atom_type_(),
	atom_charge_(0.0),
	phi_(0.0),
	theta_(0.0),
	d_(0.0),
	stub_atom1_(),
	stub_atom2_(),
	stub_atom3_()
{}

/// constructor
Adduct::Adduct(
	std::string const & adduct_name,
	std::string const & atom_name,
	std::string const & atom_type_name,
	std::string const & mm_atom_type_name,
	Real const atom_charge_in,
	Real const phi_in,
	Real const theta_in,
	Real const d_in,
	std::string const & stub_atom1_name,
	std::string const & stub_atom2_name,
	std::string const & stub_atom3_name
):
	adduct_name_( adduct_name ),
	atom_name_( atom_name ),
	atom_type_name_( atom_type_name ),
	//atom_type_(),
	mm_atom_type_name_( mm_atom_type_name ),
	//mm_atom_type_(),
	atom_charge_( atom_charge_in ),
	phi_( phi_in ),
	theta_( theta_in ),
	d_( d_in ),
	stub_atom1_( stub_atom1_name ),
	stub_atom2_( stub_atom2_name ),
	stub_atom3_( stub_atom3_name )
{}

/// accessor to adduct_name string
std::string const &
Adduct::adduct_name() const
{
	return adduct_name_;
}

/// accessor to atom_name string
std::string const &
Adduct::atom_name() const
{
	return atom_name_;
}

/// accessor to atom type string
std::string const &
Adduct::atom_type_name() const
{
	return atom_type_name_;
}

/// accessor to mm type string
std::string const &
Adduct::mm_atom_type_name() const
{
	return mm_atom_type_name_;
}

Real
Adduct::atom_charge() const
{
	return atom_charge_;
}

/// accessor for Adduct geometric info
Real
Adduct::phi() const
{
	return phi_;
}


Real
Adduct::theta() const
{
	return theta_;
}


Real
Adduct::d() const
{
	return d_;
}

/// accessor to stub_atom1 name string
std::string const &
Adduct::stub_atom1() const
{
	return stub_atom1_;
}

/// accessor to stub_atom2 name string
std::string const &
Adduct::stub_atom2() const
{
	return stub_atom2_;
}

/// accessor to stub_atom3 name string
std::string const &
Adduct::stub_atom3() const
{
	return stub_atom3_;
}

/// const accessor to stub_atom strings by index
std::string const &
Adduct::stub_atom( int const atm ) const
{
	switch( atm ) {
	case 1 : return stub_atom1_;
	case 2 : return stub_atom2_;
	case 3 : return stub_atom3_;
	}
	utility_exit_with_message( "ICoorAtomID::stub_atom should be 1--3" );
	return stub_atom1_;
}

} // chemical
} // core



#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::Adduct::save( Archive & arc ) const {
	arc( CEREAL_NVP( adduct_name_ ) ); // std::string
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
	arc( CEREAL_NVP( atom_type_name_ ) ); // std::string
	arc( CEREAL_NVP( mm_atom_type_name_ ) ); // std::string
	arc( CEREAL_NVP( atom_charge_ ) ); // Real
	arc( CEREAL_NVP( phi_ ) ); // Real
	arc( CEREAL_NVP( theta_ ) ); // Real
	arc( CEREAL_NVP( d_ ) ); // Real
	arc( CEREAL_NVP( stub_atom1_ ) ); // std::string
	arc( CEREAL_NVP( stub_atom2_ ) ); // std::string
	arc( CEREAL_NVP( stub_atom3_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::Adduct::load( Archive & arc ) {
	arc( adduct_name_ ); // std::string
	arc( atom_name_ ); // std::string
	arc( atom_type_name_ ); // std::string
	arc( mm_atom_type_name_ ); // std::string
	arc( atom_charge_ ); // Real
	arc( phi_ ); // Real
	arc( theta_ ); // Real
	arc( d_ ); // Real
	arc( stub_atom1_ ); // std::string
	arc( stub_atom2_ ); // std::string
	arc( stub_atom3_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::Adduct );
#endif // SERIALIZATION
