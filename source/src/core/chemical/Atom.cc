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
#include <core/chemical/AtomProperties.hh>

// Basic header
#include <basic/Tracer.hh>

// Utility header
#include <utility/exit.hh>

#ifdef    SERIALIZATION
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>
#include <utility/pointer/memory.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

static basic::Tracer TR( "core.chemical.Atom" );

/// @details All its properties are un-set by default.
Atom::Atom() :
	properties_( utility::pointer::make_shared< AtomProperties >() )
{}

/// @details Rosetta AtomTypes should be set through the ResidueType to ensure data consistency.
Atom::Atom(
	std::string const & name_in,
	std::string const & mm_name,
	ElementCOP element,
	Real const charge,
	Vector const & ideal_xyz ) :
	name_( name_in ),
	mm_name_( mm_name ),
	element_(std::move( element )),
	charge_( charge ),
	ideal_xyz_( ideal_xyz ),
	properties_( utility::pointer::make_shared< AtomProperties >() )
{}

bool
Atom::operator==( Atom const & atom ) const {
	return name_== atom.name_ &&
		mm_name_ == atom.mm_name_ &&
		atom_type_index_ == atom.atom_type_index_ &&
		element_ == atom.element_ &&
		gasteiger_atom_type_ == atom.gasteiger_atom_type_ &&
		is_backbone_ == atom.is_backbone_ &&
		is_actcoord_ == atom.is_actcoord_ &&
		formal_charge_ == atom.formal_charge_ &&
		charge_ == atom.charge_ &&
		abs_stereochem_ == atom.abs_stereochem_ &&
		greek_d_ == atom.greek_d_ &&
		ideal_xyz_ == atom.ideal_xyz_ &&
		( icoor_ == atom.icoor_ || ( icoor_ != nullptr && atom.icoor_ != nullptr && *icoor_ == *atom.icoor_ ) ) &&
		*properties_ == *atom.properties_ &&
		bonded_orbitals_ == atom.bonded_orbitals_;
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

void Atom::gasteiger_atom_type( core::chemical::gasteiger::GasteigerAtomTypeDataCOP gasteiger_atom_type ) {
	gasteiger_atom_type_ = gasteiger_atom_type;
}

void
Atom::set_absolute_stereochemistry( char const setting )
{
	if ( !( setting == 'R' || setting == 'S' || setting == char() ) ) {
		utility_exit_with_message( "Atom::set_absolute_stereochemistry: "
			"R and S are the only valid deginations for absolute stereochemistry." );
	}
	abs_stereochem_ = setting;
}

bool
Atom::has_property( std::string const & property ) const
{
	return properties_->has_property( property );
}

bool
Atom::has_property( AtomProperty const property ) const
{
	return properties_->has_property( property );
}

/// @brief  Generic property setting.
void
Atom::set_property( std::string const & property, bool const setting)
{
	properties_->set_property( property, setting );
}

void
Atom::set_property( AtomProperty const property, bool const setting)
{
	properties_->set_property( property, setting );
}

void
Atom::reset_all_properies(AtomProperties const & setting)
{
	properties_ = AtomPropertiesOP( new AtomProperties( setting ) );
}

// Return true if this represents a fake/mock atom.
/// @note  (Real atoms have elements that exist on the periodic table.)
bool
Atom::is_fake() const {
	if ( is_virtual() ) { return true; }
	if ( element_ ) {
		return element_->is_fake();
	} else {
		TR.Warning << "Attempted to determine real/fake status of atom without an element: " << name() << std::endl;
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
	out << "CHARMm: " << mm_name_  << "; ";
	out << "Gasteiger: ";
	if ( gasteiger_atom_type() ) {
		out << gasteiger_atom_type_->get_name();
	} else {
		out << "None";
	}
	out << endl;

	if ( is_backbone() ) {
		out << "   Backbone Atom" << endl;
	}
	out << "   Charge: " << "partial: " << charge_ << "; formal: ";
	if ( formal_charge_ > 0 ) {
		out << '+';
	}
	out << formal_charge_ << endl;
	if ( abs_stereochem_ != char() ) {
		out << "   Stereochemistry: " << abs_stereochem_ << endl;
	}
	if ( greek_d_ != NA_GREEK_DISTANCE ) {
		out << "   Greek Distance: " << greek_d_ << endl;
	}

	out << "   Properties: " << endl;
	out << *properties_;

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

#ifdef    SERIALIZATION

void
core::chemical::Atom::update_typesets( MutableResidueType const & parent ) {
	element_ = nullptr;
	if ( parent.element_set_ptr() ) {
		element_ = parent.element_set_ptr()->element( element_enum_ );
	}
	gasteiger_atom_type_ = nullptr;
	if ( ! gasteiger_atom_type_name_.empty() && parent.gasteiger_atom_typeset() ) {
		gasteiger_atom_type_ =  parent.gasteiger_atom_typeset()->atom_type( gasteiger_atom_type_name_ );
	}

}

template< class Archive >
void
core::chemical::Atom::save( Archive & arc ) const {
	arc( CEREAL_NVP( name_ ) ); // std::string
	arc( CEREAL_NVP( mm_name_ ) ); // std::string
	arc( CEREAL_NVP( atom_type_index_ ) ); // Size
	// EXEMPT element_ element_enum_
	arc( CEREAL_NVP_( "element_enum_", element() )); // core::chemical::element::Elements
	// EXEMPT gasteiger_atom_type_ gasteiger_atom_type_name_
	std::string gasteiger_atom_type_name("");
	if ( gasteiger_atom_type_ ) {
		gasteiger_atom_type_name = gasteiger_atom_type_->get_name();
	}
	arc( CEREAL_NVP_( "gasteiger_atom_type_name_", gasteiger_atom_type_name ) ); // std::string
	arc( CEREAL_NVP( is_backbone_ ) ); // bool
	arc( CEREAL_NVP( is_actcoord_ ) ); // bool
	arc( CEREAL_NVP( formal_charge_ ) ); // int
	arc( CEREAL_NVP( charge_ ) ); // Real
	arc( CEREAL_NVP( abs_stereochem_ ) ); // char
	arc( CEREAL_NVP( greek_d_ ) ); // enum core::chemical::GreekDistance
	arc( CEREAL_NVP( ideal_xyz_ ) ); // Vector
	arc( CEREAL_NVP( icoor_ ) ); // MutableICoorRecordCOP
	arc( CEREAL_NVP( bonded_orbitals_ ) ); // utility::vector1< Size >
	arc( properties_.get_op() ); // AtomPropertiesOP
}

template< class Archive >
void
core::chemical::Atom::load( Archive & arc ) {
	arc( name_ ); // std::string
	arc( mm_name_ ); // std::string
	arc( atom_type_index_ ); // Size
	// EXEMPT element_
	arc( element_enum_ ); // core::chemical::element::Elements
	// EXEMPT gasteiger_atom_type_
	arc( gasteiger_atom_type_name_ );
	arc( is_backbone_ ); // bool
	arc( is_actcoord_ ); // bool
	arc( formal_charge_ ); // int
	arc( charge_ ); // Real
	arc( abs_stereochem_ ); // char
	arc( greek_d_ ); // enum core::chemical::GreekDistance
	arc( ideal_xyz_ ); // Vector
	MutableICoorRecordOP icoor; // EXEMPT icoor_
	arc( icoor ); // MutableICoorRecordCOP
	icoor_ = icoor;
	arc( bonded_orbitals_ ); // utility::vector1< Size >
	AtomPropertiesOP props;
	arc( props ); // AtomPropertiesOP
	properties_ = props; // EXEMPT properties_
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::Atom );
#endif // SERIALIZATION
