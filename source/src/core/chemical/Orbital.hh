// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A class for defining chemical atoms, with properties specific to a ResidueType, not conformation info
/// specific to a Residue. Conformation info goes in conformaton::Orbital. OrbitalTypes are not ResidueType specific.
///
///
///
///
/// @author
/// Gordon Lemmon
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_core_chemical_Orbital_hh
#define INCLUDED_core_chemical_Orbital_hh


// Unit headers
#include <core/chemical/Orbital.fwd.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>

// Package headers
#include <core/chemical/types.hh>

// Utility headers
#include <utility/vector1_bool.hh>

// C++ headers
#include <iosfwd>

namespace core {
namespace chemical {

/// @brief basic chemical atom
///
/// @details name, element, certain properties and parameters from .params file
///
class Orbital {

public:

	/// @brief Construct a new atom type with its name and element.
	///
	/// @details All its properties are unset by default.
	///
	Orbital():
		name_(""),
		orbital_type_index_(0),
		ideal_xyz_(),
		icoor_(),
		new_icoor_()
	{}

	Orbital(
		std::string const & name_in,
		Size const orbital_type_index,
		Vector const & xyz
	):
		name_( name_in ),
		orbital_type_index_(orbital_type_index),
		ideal_xyz_(xyz),
		icoor_(),
		new_icoor_()
	{}

	Orbital(Orbital const & src) :
		name_( src.name_ ),
		orbital_type_index_(src.orbital_type_index_),
		ideal_xyz_(src.ideal_xyz_),
		icoor_(src.icoor_),
		new_icoor_(src.new_icoor_)
	{}

	void
	print( std::ostream & out ) const;

	friend
	std::ostream &
	operator<< ( std::ostream & out, const Orbital & atom_type );

	// Const Getters
	std::string const& name() const { return name_; };
	Size const& orbital_type_index() const { return orbital_type_index_; };
	Vector const& ideal_xyz() const { return ideal_xyz_; };
	orbitals::ICoorOrbitalData const& icoor() const { return icoor_; };
	orbitals::ICoorOrbitalData const& new_icoor() const { return new_icoor_; };
	// Non-const getters
	orbitals::ICoorOrbitalData & icoor() { return icoor_; };
	orbitals::ICoorOrbitalData & new_icoor() { return new_icoor_; };
	// Setters
	void name( std::string const & name ) { name_ = name; };
	void orbital_type_index( Size const & atom_type_index ) { orbital_type_index_ = atom_type_index; };
	void ideal_xyz( Vector const & xyz) { ideal_xyz_= xyz; };
	void icoor( orbitals::ICoorOrbitalData const & icoor) { icoor_ = icoor; };
	void new_icoor( orbitals::ICoorOrbitalData const & new_icoor) { new_icoor_ = new_icoor; };

	// data
private:

	// Primary data
	std::string name_;
	Size orbital_type_index_;
	// ideal xyz coordinates
	Vector ideal_xyz_;
	// ideal internal coordinates
	orbitals::ICoorOrbitalData icoor_;
	orbitals::ICoorOrbitalData new_icoor_;
};


} // chemical
} // core


#endif // INCLUDED_core_chemical_Orbital_HH
