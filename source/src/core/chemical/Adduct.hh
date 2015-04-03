// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief declaration of implementation class for residue adducts
/// @author Jim Havranek


#ifndef INCLUDED_core_chemical_Adduct_hh
#define INCLUDED_core_chemical_Adduct_hh

// Unit headers
#include <core/chemical/Adduct.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>

#include <string>

namespace core {
namespace chemical {

/// @brief Description of optional single-atom residue adducts
class Adduct {
public:
	/// default constructor
	Adduct();

	/// constructor
	Adduct(
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
	);

public:
	/// accessor to adduct_name string
	std::string const & adduct_name() const;

	/// accessor to atom_name string
	std::string const & atom_name() const;

	/// accessor to atom type string
	std::string const & atom_type_name() const;

	/// accessor to mm type string
	std::string const & mm_atom_type_name() const;

	Real atom_charge() const;

	/// @brief accessor for Adduct geometric info
	Real phi() const;


	Real theta() const;


	Real d() const;

	/// accessor to stub_atom1 name string
	std::string const & stub_atom1() const;

	/// accessor to stub_atom2 name string
	std::string const & stub_atom2() const;

	/// accessor to stub_atom3 name string
	std::string const & stub_atom3() const;

	/// const accessor to stub_atom strings by index
	std::string const & stub_atom( int const atm ) const;

private:
	std::string adduct_name_;
	std::string atom_name_;
	std::string atom_type_name_;
	// unused AtomType* atom_type_;
	std::string mm_atom_type_name_;
	// unused MMAtomType *mm_atom_type_;
	Real atom_charge_;
	Real phi_;
	Real theta_;
	Real d_;
	std::string stub_atom1_;
	std::string stub_atom2_;
	std::string stub_atom3_;
};

} // chemical
} // core


#endif
