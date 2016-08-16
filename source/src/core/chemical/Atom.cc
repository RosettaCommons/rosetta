// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/Atom.hh
/// @brief  Method definitions for chemical::Atom
/// @note   not to be confused with conformation::Atom
/// @author Phil Bradley
/// @author Labonte <JWLabonte@jhu.edu>

// Unit header
#include <core/chemical/Atom.hh>

// Package header
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>

// Basic header
#include <basic/Tracer.hh>

// Utility header
#include <utility/exit.hh>

namespace core {
namespace chemical {

static THREAD_LOCAL basic::Tracer TR( "core.chemical.Atom" );

/// @details All its properties are un-set by default.
Atom::Atom() :
	name_( "" ),
	mm_name_( "" ),
	atom_type_index_( 0 ),
	mm_atom_type_index_( 0 ),
	element_( /* 0 */ ),
	gasteiger_atom_type_( /* 0 */ ),
	formal_charge_( 0 ),
	charge_( 0 ),
	ideal_xyz_(),
	properties_( AtomPropertiesOP( new AtomProperties() ) ),
	is_hydrogen_( false ),
	has_orbitals_( false ),
	abs_stereochem_( char() ),
	greek_d_( NA_GREEK_DISTANCE )
{}

/// @details Rosetta AtomTypes should be set through the ResidueType to ensure data consistency.
Atom::Atom(
	std::string const & name_in,
	std::string const & mm_name,
	Size const mm_atom_type_index,
	ElementCOP element,
	Real const charge,
	Vector const & ideal_xyz ) :
	name_( name_in ),
	mm_name_( mm_name ),
	atom_type_index_( 0 ),
	mm_atom_type_index_( mm_atom_type_index ),
	element_( element ),
	gasteiger_atom_type_( /* 0 */ ),
	formal_charge_( 0 ),
	charge_( charge ),
	ideal_xyz_( ideal_xyz ),
	properties_( AtomPropertiesOP( new AtomProperties() ) ),
	is_hydrogen_( false ),
	has_orbitals_( false ),
	abs_stereochem_( char() ),
	greek_d_( NA_GREEK_DISTANCE )
{}

Atom::Atom( Atom const & src ) :
	name_( src.name_ ),
	mm_name_( src.mm_name_ ),
	atom_type_index_( src.atom_type_index_ ),
	mm_atom_type_index_( src.mm_atom_type_index_ ),
	element_( src.element_ ),
	gasteiger_atom_type_( src.gasteiger_atom_type_ ),
	formal_charge_( src.formal_charge_ ),
	charge_( src.charge_ ),
	ideal_xyz_( src.ideal_xyz_ ),
	properties_( AtomPropertiesOP ( new AtomProperties( *src.properties_ ) ) ),
	is_hydrogen_( 0 ),
	has_orbitals_( 0 ),
	bonded_orbitals_(),
	abs_stereochem_( src.abs_stereochem_ ),
	greek_d_( src.greek_d_ )
{}

Atom::~Atom() {}

//because you have an owning pointer in private member data, you need to implement an = operator
Atom &
Atom::operator=( Atom const & rhs )
{
	name_= rhs.name_;
	mm_name_ = rhs.mm_name_;
	atom_type_index_ = rhs.atom_type_index_;
	mm_atom_type_index_ = rhs.mm_atom_type_index_;
	element_ = rhs.element_;
	formal_charge_ = rhs.formal_charge_;
	charge_ = rhs.charge_;
	ideal_xyz_ = rhs.ideal_xyz_;
	properties_ = AtomPropertiesOP ( new AtomProperties( *rhs.properties_ ) );
	gasteiger_atom_type_ = rhs.gasteiger_atom_type_;
	is_hydrogen_ = rhs.is_hydrogen_;
	has_orbitals_ = rhs.has_orbitals_;
	bonded_orbitals_ = rhs.bonded_orbitals_;
	abs_stereochem_ = rhs.abs_stereochem_;
	greek_d_ = rhs.greek_d_;
	return *this;
}

core::chemical::element::Elements
Atom::element() const {
	if ( ! element_ ) {
		return core::chemical::element::UnknownElement;
	} else {
		return element_->element();
	}
}

gasteiger::GasteigerAtomTypeDataCOP Atom::gasteiger_atom_type() const { return gasteiger_atom_type_; }

void Atom::gasteiger_atom_type( core::chemical::gasteiger::GasteigerAtomTypeDataCOP gasteiger_atom_type ) { gasteiger_atom_type_ = gasteiger_atom_type; }

void
Atom::set_absolute_stereochemistry( char const setting )
{
	if ( !( setting == 'R' || setting == 'S' || setting == char() ) ) {
		utility_exit_with_message( "Atom::set_absolute_stereochemistry: "
			"R and S are the only valid deginations for absolute stereochemistry." );
	}
	abs_stereochem_ = setting;
}

void
Atom::greek_distance( GreekDistance const setting )
{
	greek_d_ = setting;
}

// Return true if this represents a fake/mock atom.
/// @note  (Real atoms have elements that exist on the periodic table.)
bool
Atom::is_fake() const {
	if ( is_virtual() ) { return true; }
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
	if ( ! is_virtual() ) {
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
	if ( heavyatom_has_polar_hydrogens() ) {
		out << "H-bond donor, ";
	}
	if ( is_acceptor() ) {
		out << "H-bond acceptor, ";
	}
	if ( is_hydrogen_ ) {
		if ( is_polar_hydrogen() ) {
			out << "polar hydrogen, ";
		} else if ( is_haro() ) {
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
