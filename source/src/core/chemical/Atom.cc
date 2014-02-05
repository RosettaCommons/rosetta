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
/// A class for defining atom parameters, known as atom_types
///
/// @details
/// This class contains the "chemical" information for atoms. This does not contain the actual
/// xyz coordinates of the class (xyz found in core/conformation/Atom.hh. The atom_type properties
/// are assigned by the class AtomSet which is initiated from the ChemicalManager. Atom type properties
/// are currently are read in from the file located chemical/atom_type_sets/fa_standard/atom_properties.txt.
/// These properties contain the the properties of LJ_RADIUS, LJ_WDEPTH, LK_DGRFREE, LK_LAMBDA, LK_VOLUME.
/// These properties are used in the scoring function fa_atr, fa_rep, fa_sol, which is located in the Etable (core/scoring/etable/Etable.hh)
/// Additional parameters are acceptor/donor, hybridzation, and orbital paramaters.
///
///
///
/// @author
/// Phil Bradley
/// Steven Combs - comments
///
///
/// @last_modified December 6 2010
/////////////////////////////////////////////////////////////////////////

// Rosetta headers
#include <core/chemical/Atom.hh>

#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>

// Utility headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <algorithm>

namespace core {
namespace chemical {

Atom::Atom():
		name_(""),
		mm_name_(""),
		atom_type_index_(0),
		mm_atom_type_index_(0),
		element_(0),
		gasteiger_atom_type_(0),
		formal_charge_(0),
		charge_(0),
		ideal_xyz_(),
		is_acceptor_(0),
		is_polar_hydrogen_(0),
		is_hydrogen_(0),
		is_haro_(0),
		is_virtual_(0),
		has_orbitals_(0)
{}

Atom::Atom(
		std::string const & name_in,
		std::string const mm_name,
		Size const atom_type_index,
		Size const mm_atom_type_index,
		ElementCOP element,
		Real const charge,
		Vector const ideal_xyz

):
	name_( name_in ),
	mm_name_(mm_name),
	atom_type_index_(atom_type_index),
	mm_atom_type_index_(mm_atom_type_index),
	element_(element),
	gasteiger_atom_type_(0),
	formal_charge_(0),
	charge_(charge),
	ideal_xyz_(ideal_xyz)
{}

Atom::Atom(Atom const & src) :
	name_( src.name_ ),
	mm_name_(src.mm_name_),
	atom_type_index_(src.atom_type_index_),
	mm_atom_type_index_(src.mm_atom_type_index_),
	element_(src.element_),
	gasteiger_atom_type_(src.gasteiger_atom_type_),
	formal_charge_(src.formal_charge_),
	charge_(src.charge_),
	ideal_xyz_(src.ideal_xyz_),
	heavyatom_has_polar_hydrogens_(0),
	is_acceptor_(0),
	is_polar_hydrogen_(0),
	is_hydrogen_(0),
	is_haro_(0),
	is_virtual_(0),
	has_orbitals_(0),
	bonded_orbitals_()
{}

Atom::~Atom(){}

//because you have an owning pointer in private member data, you need to implement an = operator
Atom & Atom::operator =(Atom const & rhs){
	name_= rhs.name_;
	mm_name_ = rhs.mm_name_;
	atom_type_index_ = rhs.atom_type_index_;
	mm_atom_type_index_ = rhs.mm_atom_type_index_;
	element_ = rhs.element_;
	formal_charge_ = rhs.formal_charge_;
	charge_ = rhs.charge_;
	ideal_xyz_ = rhs.ideal_xyz_;
	gasteiger_atom_type_ = rhs.gasteiger_atom_type_;
	heavyatom_has_polar_hydrogens_ = rhs.heavyatom_has_polar_hydrogens_;
	is_acceptor_ = rhs.is_acceptor_;
	is_polar_hydrogen_ = rhs.is_polar_hydrogen_;
	is_hydrogen_ = rhs.is_hydrogen_;
	is_haro_ = rhs.is_haro_;
	is_virtual_ = rhs.is_virtual_;
	has_orbitals_ = rhs.has_orbitals_;
	bonded_orbitals_ = rhs.bonded_orbitals_;
	return *this;
}

gasteiger::GasteigerAtomTypeDataCOP Atom::gasteiger_atom_type() const { return gasteiger_atom_type_; }

void Atom::gasteiger_atom_type( core::chemical::gasteiger::GasteigerAtomTypeDataCOP gasteiger_atom_type ) { gasteiger_atom_type_ = gasteiger_atom_type; }

void
Atom::print(
	std::ostream & out
) const {
	out << "Name: " << name() << std::endl;
//	out << "Type: " << type() << std::endl;
	out << "MM Name: " << mm_name() << std::endl;
	out << "atom_type_index: " << atom_type_index() << std::endl;
	out << "mm_atom_type_index: " << mm_atom_type_index() << std::endl;
	if( gasteiger_atom_type() ) {
		out << "gasteiger_atom_type: " << gasteiger_atom_type()->get_name() << std::endl;
	} else {
		out << "gasteiger_atom_type: (None)" << std::endl;
	}
	out << "charge: " << charge() << std::endl;
	out << std::endl;
}


std::ostream &
operator<< (std::ostream & out, Atom const & atom ){
	atom.print( out );
	return out;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


} // pose
} // core
