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
#include <utility/exit.hh>

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


