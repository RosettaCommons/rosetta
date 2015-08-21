// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/Atom.hh
/// @brief  Method definitions for chemical::Atom
/// @note   not to be confused with conformation::Atom
/// @author Phil Bradley

// Unit header
#include <core/chemical/Atom.hh>

// Package header
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>

// Basic header
#include <basic/Tracer.hh>

// Utility header
#include <utility/exit.hh>

// C++ header
#include <algorithm>


namespace core {
namespace chemical {

static thread_local basic::Tracer TR( "core.chemical.Atom" );

/// @details All its properties are unset by default.
Atom::Atom():
	name_(""),
	mm_name_(""),
	atom_type_index_(0),
	mm_atom_type_index_(0),
	element_(/* 0 */),
	gasteiger_atom_type_(/* 0 */),
	formal_charge_(0),
	charge_(0),
	ideal_xyz_(),
	heavyatom_has_polar_hydrogens_( false ),
	is_acceptor_( false ),
	is_polar_hydrogen_( false ),
	is_hydrogen_( false ),
	is_haro_( false ),
	is_virtual_( false ),
	has_orbitals_( false )
{}

/// @details Rosetta AtomTypes should be set through the ResidueType to ensure data consistency.
// Do NOT change to pass by reference
Atom::Atom(
	std::string const name_in,
	std::string const mm_name,
	Size const mm_atom_type_index,
	ElementCOP element,
	Real const charge,
	Vector const & ideal_xyz
):
	name_( name_in ),
	mm_name_(mm_name),
	atom_type_index_(0),
	mm_atom_type_index_(mm_atom_type_index),
	element_(element),
	gasteiger_atom_type_(/* 0 */),
	formal_charge_(0),
	charge_(charge),
	ideal_xyz_(ideal_xyz),
	heavyatom_has_polar_hydrogens_( false ),
	is_acceptor_( false ),
	is_polar_hydrogen_( false ),
	is_hydrogen_( false ),
	is_haro_( false ),
	is_virtual_( false ),
	has_orbitals_( false )
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

// Return true if this represents a fake/mock atom.
/// @note  (Real atoms have elements that exist on the periodic table.)
bool
Atom::is_fake() const {
	if ( is_virtual_ ) { return true; }
	if ( element_ ) {
		return element_->is_fake();
	} else {
		TR.Warning << "Warning: Attempted to determine real/fake status of atom without an element: " << name() << std::endl;
		return true; // Can't be real if it doesn't have an element.
	}
}

void
Atom::show( std::ostream & out ) const
{
	using namespace std;

	out << name_;
	if ( ! is_virtual_ ) {
		out << " (" << element_->get_chemical_symbol() << ')';
	} else {
		out << " (virtual)";
	}
	out << endl;

	out << "   Types (type set indices): ";
	out << "Rosetta: " /*<< type_name_*/ << " (" << atom_type_index_ << "); ";
	out << "CHARMm: " << mm_name_ << " (" << mm_atom_type_index_ << "); ";
	out << "Gasteiger: ";
	if ( gasteiger_atom_type() ) {
		out << gasteiger_atom_type_->get_name();
	} else {
		out << "None";
	}
	out << endl;

	out << "   Charge: " << "partial: " << charge_ << "; formal: ";
	if ( formal_charge_ > 0 ) {
		out << '+';
	}
	out << formal_charge_ << endl;

	out << "   Properties: ";
	if ( heavyatom_has_polar_hydrogens_ ) {
		out << "H-bond donor, ";
	}
	if ( is_acceptor_ ) {
		out << "H-bond acceptor, ";
	}
	if ( is_hydrogen_ ) {
		if ( is_polar_hydrogen_ ) {
			out << "polar hydrogen, ";
		} else if ( is_haro_ ) {
			out << "aromatic hydrogen, ";
		} else {
			out << "non-polar hydrogen, ";
		}
	}
	out << std::endl;
}


std::ostream &
operator<< ( std::ostream & out, Atom const & atom )
{
	atom.show( out );
	return out;
}

} // pose
} // core
