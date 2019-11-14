// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/Atom.hh
/// @brief  Class definitions for chemical::Atom
/// @note   not to be confused with conformation::Atom
/// @author Steven Combs

#ifndef INCLUDED_core_chemical_Atom_HH
#define INCLUDED_core_chemical_Atom_HH

// Unit headers
#include <core/chemical/Atom.fwd.hh>
#include <core/chemical/AtomProperty.hh>
#include <core/chemical/AtomProperties.fwd.hh>
#include <core/chemical/MutableICoorRecord.hh>

// Package headers
#include <core/chemical/GreekDistance.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.fwd.hh>
#include <core/chemical/Element.hh>

// Project headers
#include <core/types.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1_bool.hh>
#include <utility/pointer/deep_copy.hh>

// C++ headers
#include <string>


#ifdef    SERIALIZATION
#include <core/chemical/ResidueType.fwd.hh>
#endif

namespace core {
namespace chemical {

/// @details This class contains the "chemical" information for atoms in a MutableResidueType.
/// This does not contain the actual xyz coordinates of the atom, (which are found in core/conformation/Atom.hh)
/// It is also not used for the plain ResidueType class, which holds the corresponding information intenrally.
/// This class should contain the information that's associated with the atom,
/// calculated from other info. (Do that in the MutableResidueType -> ResidueType transition. )
/// @note    chemical::Atoms are stored in chemical::MutableResidueType (within its ResidueGraph);
/// conformation::Atoms are stored in conformation::Residue
class Atom {

public:
	/// @brief Construct a new atom type with its name and element.
	Atom();

	/// @brief Construct a new atom with the name, mm type, element, charge and position.
	Atom(
		std::string const & name_in,
		std::string const & mm_name,
		ElementCOP element,
		Real const charge,
		Vector const & ideal_xyz );

	bool operator==( Atom const & atom ) const;

	////////////////
	// Const Getters

	std::string const& name() const { return name_; }
	std::string const& mm_name() const { return mm_name_; }
	Size const& atom_type_index() const { return atom_type_index_; }
	ElementCOP element_type() const {return element_;}

	/// @brief Convenience function to go directly to the element enum
	core::chemical::element::Elements element() const;

	gasteiger::GasteigerAtomTypeDataCOP gasteiger_atom_type() const;
	bool const& is_backbone() const { return is_backbone_; }
	/// @brief Is this atom part of the action coordinate centers?
	/// @details The geometric center of all atoms listed as actcoords are the residue's "action coordinate" (for pair energy)
	bool const& is_actcoord() const { return is_actcoord_; }
	int const& formal_charge() const { return formal_charge_; }
	Real const& charge() const { return charge_; }
	Vector const& ideal_xyz() const { return ideal_xyz_; }

	/// @brief The ICoor record for this residue -- may be null.
	MutableICoorRecordCOP icoor() const { return icoor_; }

	/// @brief  Return the absolute stereochemistry (R/S designation) of this stereocenter.
	char absolute_stereochemistry() const { return abs_stereochem_; }

	/// @brief  How far (in Greek letters) is this atom from the primary functional group of the molecule?
	GreekDistance greek_distance() const { return greek_d_; }

	/// @brief Which orbital indicies are attached to this atom?
	utility::vector1< core::Size > const &
	bonded_orbitals() const { return bonded_orbitals_; }

	//////////
	// Setters

	void element_type(ElementCOP element) {element_ = element;}
	void mm_name( std::string const & name ) { mm_name_ = name; };
	void gasteiger_atom_type( core::chemical::gasteiger::GasteigerAtomTypeDataCOP gasteiger_atom_type );

	void is_backbone( bool setting ) { is_backbone_ = setting; }
	void is_actcoord( bool setting ) { is_actcoord_ = setting; }
	void formal_charge( int charge ) { formal_charge_ = charge; }
	void charge( Real const & charge ) { charge_ = charge; };
	void ideal_xyz( Vector const & ideal_xyz ) {
		ideal_xyz_= ideal_xyz;
	}
	void icoor( MutableICoorRecordCOP icoor_record ) { icoor_ = icoor_record; }

	/// @brief  Set the absolute stereochemistry (R/S designation) of this stereocenter.
	void set_absolute_stereochemistry( char const setting );

	/// @brief  Set how far (in Greek letters) this atom is from the primary functional group of the molecule.
	void greek_distance( GreekDistance const setting ) { greek_d_ = setting; }

	void set_bonded_orbitals( utility::vector1<core::Size> const & setting ) {
		bonded_orbitals_ = setting;
	}

	void add_bonded_orbital( Size orbital_index ) {
		bonded_orbitals_.push_back( orbital_index );
	}

	/////////////////////
	// Protected Setters
	//
	// These are protected because they should really only be changes through the owning MutableResidueType
	// (As there's certain bookkeeping that the MutableResidueType needs to do.)
protected:
	friend MutableResidueType; // Only for these functions! don't access other private members!

	// @details Use MutableResidueType::rename_atom() instead.
	void name( std::string const & name ) { name_ = name; };
	// @details Use MutableResidueType::set_atom_type() instead.
	void atom_type_index( Size const & atom_type_index ) { atom_type_index_ = atom_type_index; };

public:
	//////////////
	// Properties

	/// @brief Access the collection of properties for this Atom.
	/// IMPORTANT -- This only looks at manually set properties, not any automatically derived ones.
	AtomProperties const & properties() const { return *properties_; }

	/// @brief  Generic property access.
	/// IMPORTANT -- This only looks at manually set properties, not any automatically derived ones.
	bool
	has_property( std::string const & property ) const;

	/// @brief  Generic property access.
	/// IMPORTANT -- This only looks at manually set properties, not any automatically derived ones.
	bool
	has_property( AtomProperty const property ) const;

	//teiger_atom_type @brief  Generic property setting.
	void
	set_property( std::string const & property, bool const setting);

	void
	set_property( AtomProperty const property, bool const setting);

	/// @brief Discard all the current properties, and set the the passed values.
	void
	reset_all_properies(AtomProperties const & setting);

	/// @brief  Set whether or not this atom is virtual.
	void is_virtual( bool setting ) { set_property( VIRTUAL_ATOM, setting ); }

	////////////////////
	// Derived info

	bool is_hydrogen() const { return element() == element::H; }

	/// @brief  Get whether or not this atom is virtual.
	bool is_virtual() const { return has_property( VIRTUAL_ATOM ); }

	/// @brief Return true if this represents a fake/mock atom.
	bool is_fake() const;

	///////////////////
	// Other functions.

	/// @brief  Generate string representation of chemical::Atom for debugging purposes.
	void show( std::ostream & out=std::cout ) const;

	// data
private:

	std::string name_;
	std::string mm_name_;
	Size atom_type_index_ = 0;
	ElementCOP element_;
	/// Gasteiger atom-type
	gasteiger::GasteigerAtomTypeDataCOP gasteiger_atom_type_;
	bool is_backbone_ = false;
	bool is_actcoord_ = false;
	/// The formal (integral) charge on the atom.
	int formal_charge_ = 0;
	Real charge_ = 0.0;
	char abs_stereochem_ = '\0';
	GreekDistance greek_d_ = NA_GREEK_DISTANCE;

	Vector ideal_xyz_{0,0,0};
	/// @brief Internal coordinates on how to build the given atom - used in preference to ideal_xyz if set.
	MutableICoorRecordCOP icoor_; // COP because it's a copy-on-write sematics.

	// General properties of this atom.
	// For the MutableResidueType Atom, only explicitly set properties are present.
	// There are other properties which are indirect consequences of the atom type which don't appear here.
	utility::pointer::DeepCopyOP< AtomProperties > properties_;

	utility::vector1< Size > bonded_orbitals_; /// @brief The orbitals bonded to this atom

#ifdef    SERIALIZATION
public:
	// WARNING: The serialization of an Atom is not 100% self-contained.
	// You need to call update_typesets on it after with the parent Restype in order for it to work.
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );

	void update_typesets( MutableResidueType const & parent );
private:
	// Temporary entries for serialization loading
	std::string gasteiger_atom_type_name_;
	core::chemical::element::Elements element_enum_;

#endif // SERIALIZATION

};

std::ostream & operator<< ( std::ostream & out, Atom const & atom);

}  // chemical
}  // core

#endif  // INCLUDED_core_chemical_Atom_HH
