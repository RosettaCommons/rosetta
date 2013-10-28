// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin Atom
///
/// @brief
/// A class for defining chemical atoms, with properties specific to a ResidueType, not conformation info
/// specific to a Residue. Conformation info goes in conformation::Atom. AtomTypes are not ResidueType specific.
///
///
///
///
/// @author
/// Gordon Lemmon
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_core_chemical_Atom_hh
#define INCLUDED_core_chemical_Atom_hh


// Unit headers
#include <core/chemical/Atom.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <core/chemical/Bond.fwd.hh> // only for Temp BondName


// Package headers
#include <core/chemical/types.hh>

// Utility headers
#include <utility/vector1_bool.hh>

// C++ headers
#include <string>

namespace core {
namespace chemical {

/// @brief basic chemical atom
///
/// @details name, element, certain properties and parameters from .params file
///
class Atom {

public:

	/// @brief Construct a new atom type with its name and element.
	///
	/// @details All its properties are unset by default.
	///
	Atom():
			name_(""),
			mm_name_(""),
			bonded_neighbors_(),
			bonded_neighbor_types_(),
			atom_type_index_(0),
			mm_atom_type_index_(0),
			charge_(0),
			ideal_xyz_(),
			icoor_()
	{}

	Atom(
			std::string const & name_in,
		//	std::string const type_name,
			std::string const mm_name,
			Size const atom_type_index,
			Size const mm_atom_type_index,
			Real const charge,
			Vector const ideal_xyz,
			AtomICoor const icoor = AtomICoor()

	):
		name_( name_in ),
		//type_name_(type_name),
		mm_name_(mm_name),
		bonded_neighbors_(),
		bonded_neighbor_types_(),
		atom_type_index_(atom_type_index),
		mm_atom_type_index_(mm_atom_type_index),
		charge_(charge),
		ideal_xyz_(ideal_xyz),
		icoor_(icoor)
	{}

	Atom(Atom const & src) :
		name_( src.name_ ),
		//type_name_(src.type_name),
		mm_name_(src.mm_name_),
		bonded_neighbors_(src.bonded_neighbors_),
		bonded_neighbor_types_(src.bonded_neighbor_types_),
		atom_type_index_(src.atom_type_index_),
		mm_atom_type_index_(src.mm_atom_type_index_),
		charge_(src.charge_),
		ideal_xyz_(src.ideal_xyz_),
		icoor_(src.icoor_)
	{}

	void
	print( std::ostream & out ) const;

	friend
	std::ostream &
	operator<< ( std::ostream & out, Atom const & atom);

	bool operator==(Atom const & atom) const{
		return	name_== atom.name_ &&
				mm_name_ == atom.mm_name_ &&
				atom_type_index_ == atom.atom_type_index_ &&
				mm_atom_type_index_ == atom.mm_atom_type_index_ &&
				charge_ == atom.charge_ &&
				bonded_neighbors_ == atom.bonded_neighbors_ &&
				bonded_neighbor_types_ == atom.bonded_neighbor_types_ &&
				ideal_xyz_ == atom.ideal_xyz_;
	}

// Const Getters
	std::string const& name() const { return name_; };
	//std::string const& type_name() const { return type_name_; };
	std::string const& mm_name() const { return mm_name_; };
	Size const& atom_type_index() const { return atom_type_index_; };
	Size const& mm_atom_type_index() const { return mm_atom_type_index_; };
	Real const& charge() const { return charge_; };
	Vector const& ideal_xyz() const { return ideal_xyz_; };
	AtomICoor const& icoor() const { return icoor_; };
    AtomIndices const& bonded_neighbors() const{ return bonded_neighbors_;}
    utility::vector1<BondName> const& bonded_neighbor_types() const{ return bonded_neighbor_types_;}
// Non-const getters
	AtomICoor & icoor() { return icoor_; };
	AtomIndices & bonded_neighbors(){ return bonded_neighbors_;}
	utility::vector1<BondName> & bonded_neighbor_types() { return bonded_neighbor_types_;}

// Setters
	void name( std::string const & name ) { name_ = name; };
	//std::string const& type_name() const { return type_name_; };
	void mm_name( std::string const & name ) { mm_name_ = name; };
	void atom_type_index( Size const & atom_type_index ) { atom_type_index_ = atom_type_index; };
	void mm_atom_type_index( Size const & mm_atom_type_index ) { mm_atom_type_index_ = mm_atom_type_index; };
	void charge( Real const & charge ) { charge_ = charge; };
	void ideal_xyz( Vector const & ideal_xyz) { ideal_xyz_= ideal_xyz; };
	void bonded_neighbors( AtomIndices const & bonded_neighbors) { bonded_neighbors_ = bonded_neighbors; }
	void icoor( AtomICoor const & icoor) { icoor_ = icoor; };

	// data
private:
	// Primary data
	std::string name_;
	//std::string const type_name_;

	// Secondary data
	std::string mm_name_;
	//	std::string const csd_atom_name_;
	AtomIndices bonded_neighbors_;
	utility::vector1<BondName> bonded_neighbor_types_;
	Size atom_type_index_;
	/// MM atom-type index
	Size mm_atom_type_index_;
	Real charge_;
	Vector ideal_xyz_;
	AtomICoor icoor_;
};


} // chemical
} // core



#endif // INCLUDED_core_chemical_Atom_HH
