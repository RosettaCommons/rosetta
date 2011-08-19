// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/Element.hh
/// @brief  Molecular mechanics atom type class
/// @author P. Douglas Renfrew (renfrew@unc.edu)


#ifndef INCLUDED_core_chemical_Element_hh
#define INCLUDED_core_chemical_Element_hh

// Unit headers
#include <core/chemical/Element.fwd.hh>

// Project headers
#include <core/types.hh>

// C++ headers
#include <string>

namespace core {
namespace chemical {

/// @brief class to describe Elements
///
/// @details class to describe elements
/// Borrows heavily and functions similarly to the rosetta atom type class, AtomType
///
class Element
{

public:

	///  @brief Construct a new Element with its name
	Element(
			Size const z,
			std::string const & symbol,
			std::string const & name,
			Real const weight,
			Size const mass
	):
		z_(z),
		symbol_(symbol),
		name_( name ),
		weight_( weight ),
		mass_( mass )
	{}
	/// @brief Return the atomic number
	core::Size z() const { return z_; }
	/// @brief Return the element symbol
	std::string const& symbol() const { return symbol_; }
	/// @brief Return the full name of the Element
	std::string const& name() const { return name_; }
	/// @brief Return the LJ radius of the atom type
	Real weight() const { return weight_; }
	/// @brief Return the LJ well depth of the atom type
	core::Size mass() const { return mass_; }

private:
	/// @brief #of protons+neutrons for most common isotope
	Size z_;

	/// @brief symbol of the element
	std::string const symbol_;

	/// @brief name of the element
	std::string const name_;

	/// @brief atomic weight (average weight of all isotopes)
	Real weight_;
	/// @brief #of protons+neutrons for most common isotope
	Size mass_;
};


} // chemical
} // core



#endif // INCLUDED_core_chemical_Element_HH
