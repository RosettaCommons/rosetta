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
/// A class for defining chemical atoms, with properties specific to the atom. This class should not have information
/// associated with the conformation or residuetype. Everything should be concerning the atom. Conformation goes
/// in core::conformation while data for ResidueTYpe is cached there
///
///
///
///
/// @author
/// Steven Combs
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_core_chemical_Atom_hh
#define INCLUDED_core_chemical_Atom_hh


// Unit headers
#include <core/chemical/Atom.fwd.hh>
//#include <core/chemical/AtomICoor.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
//#include <core/chemical/Bond.fwd.hh> // only for Temp BondName
#include <core/chemical/gasteiger/GasteigerAtomTypeData.fwd.hh>
#include <core/chemical/Element.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

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
	Atom();


	/// @brief Construct a new atom with the type it is, the mm type, and the name.
	Atom(
			std::string const & name_in,
			std::string const mm_name,
			Size const atom_type_index,
			Size const mm_atom_type_index,
			ElementCOP element,
			Real const charge,
			Vector const ideal_xyz
	);


	Atom(Atom const & src);

	//destructor (default)
	~Atom();

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
				element_ == atom.element_ &&
				formal_charge_ == atom.formal_charge_ &&
				charge_ == atom.charge_ &&
				ideal_xyz_ == atom.ideal_xyz_ &&
				gasteiger_atom_type_ == atom.gasteiger_atom_type_ &&
				heavyatom_has_polar_hydrogens_ == atom.heavyatom_has_polar_hydrogens_ &&
				is_acceptor_ == atom.is_acceptor_ &&
				is_polar_hydrogen_ == atom.is_polar_hydrogen_ &&
				is_hydrogen_ == atom.is_hydrogen_ &&
				is_haro_ == atom.is_haro_ &&
				is_virtual_ == atom.is_virtual_ &&
				has_orbitals_ == atom.has_orbitals_  &&
				bonded_orbitals_ == atom.bonded_orbitals_;
	}

	Atom & operator =(Atom const & atom);

	// Const Getters
	std::string const& name() const { return name_; }
	//std::string const& type_name() const { return type_name_; };
	std::string const& mm_name() const { return mm_name_; }
	Size const& atom_type_index() const { return atom_type_index_; }
	Size const& mm_atom_type_index() const { return mm_atom_type_index_; }
	ElementCOP element_type() const {return element_;}
	gasteiger::GasteigerAtomTypeDataCOP gasteiger_atom_type() const;
	int const& formal_charge() const { return formal_charge_; }
	Real const& charge() const { return charge_; }
	Vector const& ideal_xyz() const { return ideal_xyz_; };
	utility::vector1<Size> const & bonded_orbitals() const{return bonded_orbitals_;}
	utility::vector1<Size>  & bonded_orbitals() {return bonded_orbitals_;}
	bool heavyatom_has_polar_hydrogens() const{ return heavyatom_has_polar_hydrogens_;}
	bool is_acceptor() const{ return is_acceptor_;}
	bool is_polar_hydrogen() const{ return is_polar_hydrogen_;}
	bool is_hydrogen() const{return is_hydrogen_;}
	bool is_haro() const{return is_haro_;}
	bool is_virtual() const{return is_virtual_;}
	bool has_orbitals() const{return has_orbitals_;}
	// Non-const getters


	// Setters
	void name( std::string const & name ) { name_ = name; };
	void mm_name( std::string const & name ) { mm_name_ = name; };
	void atom_type_index( Size const & atom_type_index ) { atom_type_index_ = atom_type_index; };
	void mm_atom_type_index( Size const & mm_atom_type_index ) { mm_atom_type_index_ = mm_atom_type_index; };
	void gasteiger_atom_type( core::chemical::gasteiger::GasteigerAtomTypeDataCOP gasteiger_atom_type );
	void formal_charge( int charge ) { formal_charge_ = charge; }
	void charge( Real const & charge ) { charge_ = charge; };
	void ideal_xyz( Vector const & ideal_xyz) { ideal_xyz_= ideal_xyz; };
	void heavyatom_has_polar_hydrogens( bool heavyatom_has_polar_hydrogens){ heavyatom_has_polar_hydrogens_ = heavyatom_has_polar_hydrogens;}
	void is_polar_hydrogen(bool polar){is_polar_hydrogen_ = polar;}
	void is_hydrogen(bool hydrogen){is_hydrogen_= hydrogen;}
	void is_haro(bool haro){is_haro_ = haro;}
	void is_acceptor(bool acceptor){is_acceptor_ = acceptor;}
	void is_virtual(bool is_virtual){is_virtual_ = is_virtual;}
	void has_orbitals(bool orbitals){has_orbitals_ = orbitals;}

	// Calculated data

	/// @brief Return true if this represents a fake/mock atom.
	/// (Real atoms have elements that exist on the periodic table.)
	bool is_fake() const;

	// data
private:
	std::string name_;
	std::string mm_name_;
	Size atom_type_index_;
	Size mm_atom_type_index_;
	ElementCOP element_;
	/// Gasteiger atom-type
	gasteiger::GasteigerAtomTypeDataCOP gasteiger_atom_type_;
	/// The formal (integral) charge on the atom.
	int formal_charge_;
	Real charge_;
	Vector ideal_xyz_;
	bool heavyatom_has_polar_hydrogens_; // is an atom both a heavy atom and chemically bonded to a polar hydrogen?
	bool is_acceptor_; // is an atom both a heavy atom and capable of accepting hydrogen bonds?
	bool is_polar_hydrogen_; // is an atom a polar hydrogen?
	bool is_hydrogen_;
	bool is_haro_;
	bool is_virtual_;
	bool has_orbitals_;
	utility::vector1<Size> bonded_orbitals_;


};


} // chemical
} // core



#endif // INCLUDED_core_chemical_Atom_HH
