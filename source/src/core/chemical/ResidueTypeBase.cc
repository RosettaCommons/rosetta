// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ResidueTypeBase.cc
/// @brief Method definitions for ResidueTypeBase
/// @author
/// Phil Bradley
/// Rocco Moretti (rmorettiase@gmail.com)
/// Steven Combs
/// Vikram K. Mulligan - properties for D-, beta- and other noncanonicals
/// Jason W. Labonte (code related to rings, properties, lipids, carbohydrates, and other non-AAs)

// Unit headers
#include <core/chemical/ResidueTypeBase.hh>
#include <core/chemical/ResidueConnection.hh>

// Package Headers
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ResidueType.hh>

// Project Headers
#include <core/chemical/residue_support.hh>
#include <core/chemical/icoor_support.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/bond_support.hh>
#include <core/chemical/RestypeDestructionEvent.hh>
#include <core/chemical/rotamers/NCAARotamerLibrarySpecification.hh>
#include <core/chemical/Orbital.hh>
// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility>
#include <utility/py/PyAssert.hh>
#include <utility/vector1.hh>
#include <utility/graph/ring_detection.hh>

// External headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>

#include <numeric/conversions.hh>

// C++ headers
#include <algorithm>

#ifdef    SERIALIZATION
#include <core/chemical/ResidueGraphTypes.srlz.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/keys/Key2Tuple.srlz.hh>
#include <utility/keys/Key3Tuple.srlz.hh>
#include <utility/keys/Key4Tuple.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

static basic::Tracer TR( "core.chemical.ResidueTypeBase" );

ResidueTypeBase::ResidueTypeBase() = default; // private, not deleted because of serialization

ResidueTypeBase::ResidueTypeBase(
	AtomTypeSetCOP atom_types,
	ElementSetCOP elements,
	MMAtomTypeSetCOP mm_atom_types,
	orbitals::OrbitalTypeSetCOP orbital_types
) : utility::VirtualBase(),
	atom_types_( std::move( atom_types ) ),
	elements_( std::move( elements ) ),
	mm_atom_types_( std::move( mm_atom_types ) ),
	gasteiger_atom_types_(),
	orbital_types_( std::move( orbital_types ) ),
	properties_( utility::pointer::make_shared< ResidueProperties >() )
{
	if ( atom_types_ != nullptr ) {
		// For any actual ResidueTypeBase atom_types should be valid, but there's tricky RTS bootstraping logic
		// in Patch.cc that uses nullptr AtomTypeSets.
		mode_ = atom_types_->mode();
	}
}

//////////////////////////////////////////////////////////////////////////////

void
ResidueTypeBase::atom_type_set( AtomTypeSetCOP setting ) {
	debug_assert( setting );
	atom_types_ = setting;
	mode_ = atom_types_->mode();
}

gasteiger::GasteigerAtomTypeSetCOP
ResidueTypeBase::gasteiger_atom_typeset() const {
	return gasteiger_atom_types_;
}

void
ResidueTypeBase::set_gasteiger_atom_typeset( gasteiger::GasteigerAtomTypeSetCOP setting ) {
	gasteiger_atom_types_ = setting;
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
////////////////          Orbital Functions     //////////////////////
//////////////////////////////////////////////////////////////////////

Size
ResidueTypeBase::n_orbitals() const {
	return orbitals_.size();
}

Orbital const &
ResidueTypeBase::orbital(Size const orbital_index) const{
	return orbitals_[orbital_index];
}

Orbital const &
ResidueTypeBase::orbital(std::string const & orbital_name) const{
	return orbitals_[ orbital_index(orbital_name) ];
}

orbitals::OrbitalType const &
ResidueTypeBase::orbital_type(Size const orbital_index) const
{
	orbitals::OrbitalTypeSetCOP orbital_types( orbital_types_ );
	return ( *orbital_types )[ orbitals_[ orbital_index ].orbital_type_index() ];
}

core::Size
ResidueTypeBase::orbital_index( std::string const & name ) const
{
	auto iter( orbitals_index_.find( name ) );
	if ( iter == orbitals_index_.end() ) {
		utility_exit_with_message("unknown orbital_name: " + name3() + "  " + name );
	}
	return iter->second;
}

///////////////////////////////////////////////////////////////////////////////
// Property-Related Methods
///////////////////////////////////////////////////////////////////////////////

void
ResidueTypeBase::set_properties( ResiduePropertiesOP properties ) {
	properties_ = properties;
}

void
ResidueTypeBase::add_property( std::string const & property )
{
	core::chemical::ResidueProperty const prop_enum( core::chemical::ResidueProperties::get_property_from_string( property ) );
	runtime_assert_string_msg( prop_enum != NO_PROPERTY, "Error in ResidueType::add_property(): Could not parse property \"" + property + "\" as a valid residue type property." );
	add_property( prop_enum );
}

/// @brief Add a property to this ResidueType, by properties enum.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
ResidueTypeBase::add_property(
	core::chemical::ResidueProperty const property
) {
	runtime_assert( property != NO_PROPERTY );
	properties_->set_property( property, true );

	// AMW: as of 8/18/15, L_AA and D_AA no longer imply ALPHA_AA.

	// Special "umbrella cases"
	// FIXME: There really shouldn't be as many umbrella cases, IMO. ~Labonte
	if ( property == PROTEIN ) {
		properties_->set_property( POLYMER, true );
	} else if ( property == ALPHA_AA ) {
		properties_->set_property( PROTEIN, true );
		properties_->set_property( POLYMER, true );
	} else if ( property == BETA_AA ) {
		properties_->set_property( PROTEIN, true );
		properties_->set_property( POLYMER, true );
	} else if ( property == L_AA ) {
		properties_->set_property( PROTEIN, true );
		properties_->set_property( POLYMER, true );
		//properties_->set_property( ALPHA_AA, true );
	} else if ( property == D_AA ) {
		properties_->set_property( PROTEIN, true );
		properties_->set_property( POLYMER, true );
		//properties_->set_property( ALPHA_AA, true );
	} else if ( property == DNA ) {
		properties_->set_property( POLYMER, true );
	} else if ( property == RNA ) {
		properties_->set_property( POLYMER, true );
	} else if ( property == PEPTOID ) {
		properties_->set_property( POLYMER, true );
		// amw: This will no longer be true if we incorporate
		// alanine peptoids (or update our chirality model entirely)
		properties_->set_property( ACHIRAL_BACKBONE, true );
	} else if ( property == LOWERTERM_TRUNC ) {
		properties_->set_property( LOWER_TERMINUS, true );
	} else if ( property == UPPERTERM_TRUNC ) {
		properties_->set_property( UPPER_TERMINUS, true );
	} else if ( property == PHOSPHONATE ) {
		properties_->set_property( POLYMER, true );
	} else if ( property == PHOSPHONATE_UPPER ) {
		properties_->set_property( UPPER_TERMINUS, true );
		properties_->set_property( PHOSPHONATE, true );
	} else if ( property == ACETYLATED_NTERMINUS ) {
		properties_->set_property( LOWER_TERMINUS, true );
	} else if ( property == METHYLATED_CTERMINUS ) {
		properties_->set_property( UPPER_TERMINUS, true );
	}
}

void
ResidueTypeBase::set_adduct_flag( bool adduct_in )
{
	properties_->set_property( ADDUCT, adduct_in );
}

void
ResidueTypeBase::add_numeric_property(std::string const & tag, core::Real value)
{
	properties_->add_numeric_property( tag, value );
}

void
ResidueTypeBase::add_string_property(std::string const & tag, std::string value)
{
	properties_->add_string_property( tag, value );
}

/// @details This is needed for deleting properties, which occurs in certain PTMs.
void
ResidueTypeBase::delete_property( std::string const & property )
{
	properties_->set_property( property, false );
}

/// @brief Delete a property of this ResidueType, by properties enum.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
ResidueTypeBase::delete_property(
	core::chemical::ResidueProperty const property
) {
	runtime_assert( property != NO_PROPERTY );
	properties_->set_property( property, false );
}

/// @brief  Generic property access.
bool
ResidueTypeBase::has_property( std::string const & property ) const
{
	return properties_->has_property( property );
}

/// @brief  Generic property access, by ResidueProperty.
///
bool
ResidueTypeBase::has_property( ResidueProperty const property ) const
{
	return properties_->has_property( property );
}

core::Real
ResidueTypeBase::get_numeric_property(std::string const & tag) const
{
	std::map<std::string, core::Real> const & numeric_properties( properties_->numeric_properties() );
	auto property_it( numeric_properties.find( tag ) );
	if ( property_it == numeric_properties.end() ) {
		throw CREATE_EXCEPTION(utility::excn::KeyError,  tag + " does not exist in ResidueTypeBase with name " + name3_ );
		return 0.0; //keep compilers happy
	}

	return property_it->second;
}

std::string
ResidueTypeBase::get_string_property(std::string const & tag) const
{
	std::map<std::string, std::string> const & string_properties( properties_->string_properties() );
	auto property_it(string_properties.find(tag));
	if ( property_it == string_properties.end() ) {
		throw CREATE_EXCEPTION(utility::excn::KeyError, tag + " does not exist in ResidueTypeBase with name " + name3_);
		return "";
	}
	return property_it->second;
}


void
ResidueTypeBase::add_variant_type( VariantType const variant_type )
{
	properties_->set_variant_type( variant_type, true );
}

void
ResidueTypeBase::add_variant_type( std::string const & variant_type )
{
	properties_->set_variant_type( variant_type, true );
}

void
ResidueTypeBase::remove_variant_type( VariantType const variant_type )
{
	properties_->set_variant_type( variant_type, false );
}

void
ResidueTypeBase::remove_variant_type( std::string const & variant_type )
{
	properties_->set_variant_type( variant_type, false );
}

/// @brief  Generic variant access.
bool
ResidueTypeBase::has_variant_type( VariantType const variant_type ) const
{
	return properties_->is_variant_type( variant_type );
}

// TODO: Find a way to remove this; it only exists because of how ResidueTypeBaseSelectors are currently written. ~Labonte
/// @brief  Generic variant access by string.
bool
ResidueTypeBase::has_variant_type( std::string const & variant_type ) const
{
	return properties_->is_variant_type( variant_type );
}

/// @details "Custom" VariantTypes as strings are permitted for the enzdes and metalloproteins cases.
/// Do not enable unless you have a good reason to, as string look-ups are less efficient and more error-prone.
void
ResidueTypeBase::enable_custom_variant_types()
{
	properties_->enable_custom_variant_types();
}

utility::vector1< std::string >
ResidueTypeBase::variant_types() const {
	return properties_->get_list_of_variants();
}

/// @brief Get a vector of VariantType enums for this ResidueTypeBase.
/// @details This ONLY includes standard, enum-based variants, not on-the-fly custom variants.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
utility::vector1< VariantType >
ResidueTypeBase::variant_type_enums() const {
	return properties_->get_list_of_variant_enums();
}

/// @brief Get a list of custom VariantType strings for this ResidueTypeBase (by const reference).
/// @details This ONLY includes custom, on-the-fly variants, not standard variants.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
utility::vector1< std::string > const &
ResidueTypeBase::custom_variant_types() const{
	return properties_->get_list_of_custom_variants_by_reference();
}

////////////////////////////////////////////////////////////
///////////// Rotamer library //////////////////////////////
////////////////////////////////////////////////////////////

void
ResidueTypeBase::rotamer_library_specification( rotamers::RotamerLibrarySpecificationOP rotlibspec) {
	rotamer_library_specification_ = rotlibspec;
}

rotamers::RotamerLibrarySpecificationCOP
ResidueTypeBase::rotamer_library_specification() const {
	return rotamer_library_specification_;
}

/// @brief Nonconst access to the RotamerLibrarySpecification.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
rotamers::RotamerLibrarySpecificationOP
ResidueTypeBase::rotamer_library_specification_nonconst() {
	return rotamer_library_specification_;
}

/// @brief Remove any rotamer library specifications attached to this ResidueTypeBase.
/// @details After this operation, the rotamer_library_specification() method returns a NULL pointer.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
ResidueTypeBase::strip_rotamer_library_specification() {
	rotamer_library_specification_ = nullptr;
}

////////////////////////////////////////////////////////////////////////////
// rama methods
////////////////////////////////////////////////////////////////////////////

/// @brief Get the key name for the mainchain torsion potential.
/// @details Stored internally as a string for base residue types.  Empty string is stored by default for derived
/// residue types (in which case this function returns the string stored in the base ResidueTypeBase), though this can be overridden.
std::string const &
ResidueTypeBase::get_rama_prepro_mainchain_torsion_potential_name( bool const pre_proline_position ) const {

	if ( pre_proline_position ) {
		if ( !rama_prepro_mainchain_torsion_potential_name_beforeproline_.empty() ) return rama_prepro_mainchain_torsion_potential_name_beforeproline_;
		//If the string is empty...:
		if ( !is_base_type() ) {
			std::string const & basestring( get_base_type_cop()->get_rama_prepro_mainchain_torsion_potential_name(pre_proline_position) );
			if ( !basestring.empty() ) return basestring;
		}
		return rama_prepro_mainchain_torsion_potential_name_beforeproline_; //Returns an empty string if this is empty AND the base type is empty.
	}

	//Otherwise...

	if ( !rama_prepro_mainchain_torsion_potential_name_.empty() ) return rama_prepro_mainchain_torsion_potential_name_;
	//If the string is empty...:
	if ( !is_base_type() ) {
		std::string const & basestring( get_base_type_cop()->get_rama_prepro_mainchain_torsion_potential_name(pre_proline_position) );
		if ( !basestring.empty() ) return basestring;
	}
	return rama_prepro_mainchain_torsion_potential_name_; //Returns an empty string if this is empty AND the base type is empty.
}

/// @brief Do the rama_prepro mainchain torsion potentials of this residue match another?
///
bool
ResidueTypeBase::mainchain_potentials_match(
	ResidueTypeBase const &other
) const {
	return
		rama_prepro_mainchain_torsion_potential_name_ == other.rama_prepro_mainchain_torsion_potential_name_ &&
		rama_prepro_mainchain_torsion_potential_name_beforeproline_ == other.rama_prepro_mainchain_torsion_potential_name_beforeproline_ &&
		rama_prepro_map_file_name_ == other.rama_prepro_map_file_name_ &&
		rama_prepro_map_file_name_beforeproline_ == other.rama_prepro_map_file_name_beforeproline_
		;
}

/// @brief Set the key name for the mainchain torsion potential.
/// @details Stored internally as a string for base residue types.  Empty string is stored by default for derived
/// residue types (pointing the function to the base type), though this can be overridden using this function.
void
ResidueTypeBase::set_rama_prepro_mainchain_torsion_potential_name(
	std::string const &name_in,
	bool const pre_proline_position
) {
	if ( pre_proline_position ) {
		rama_prepro_mainchain_torsion_potential_name_beforeproline_ = name_in;
	} else {
		rama_prepro_mainchain_torsion_potential_name_ = name_in;
	}
}

/// @brief Get the file name for the mainchain torsion potential used by the RamaPrePro score term.
/// @details Stored internally as a string for base residue types.  Empty string is stored by default for derived
/// residue types (in which case this function returns the string stored in the base ResidueTypeBase), though this can be overridden.
std::string const &
ResidueTypeBase::get_rama_prepro_map_file_name( bool const pre_proline_position ) const {
	if ( pre_proline_position ) {
		if ( !rama_prepro_map_file_name_beforeproline_.empty() ) return rama_prepro_map_file_name_beforeproline_;
		//If the string is empty...:
		if ( !is_base_type() ) {
			std::string const & basestring( get_base_type_cop()->get_rama_prepro_map_file_name(pre_proline_position) );
			if ( !basestring.empty() ) return basestring;
		}
		return rama_prepro_map_file_name_beforeproline_; //Returns an empty string if this is empty AND the base type is empty.
	}

	//Otherwise...

	if ( !rama_prepro_map_file_name_.empty() ) return rama_prepro_map_file_name_;
	//If the string is empty...:
	if ( !is_base_type() ) {
		std::string const & basestring( get_base_type_cop()->get_rama_prepro_map_file_name(pre_proline_position) );
		if ( !basestring.empty() ) return basestring;
	}
	return rama_prepro_map_file_name_; //Returns an empty string if this is empty AND the base type is empty.
}

/// @brief Set the file name for the mainchain torsion potential used by the RamaPrePro score term.
/// @details Stored internally as a string for base residue types.  Empty string is stored by default for derived
/// residue types (pointing the function to the base type), though this can be overridden using this function.
void
ResidueTypeBase::set_rama_prepro_map_file_name(
	std::string const &filename_in,
	bool const pre_proline_position
) {
	if ( pre_proline_position ) {
		rama_prepro_map_file_name_beforeproline_ = filename_in;
	} else {
		rama_prepro_map_file_name_ = filename_in;
	}
}

/// @brief Returns true if and only if (a) this is not a base type, AND (b) there is a rama_prepro_mainchain_torsion_map_file_name_
/// defined for this ResidueTypeBase (which is presumably different from that of the base type).
/// @details If pre_proline_position is true, checks rama_prepro_mainchain_torsion_map_file_name_beforeproline_ instead of
/// rama_prepro_mainchain_torsion_potential_name_.
bool
ResidueTypeBase::defines_custom_rama_prepro_map( bool const pre_proline_position ) const {
	if ( is_base_type() ) return false;
	if ( pre_proline_position ) {
		return ( !rama_prepro_map_file_name_beforeproline_.empty() );
	}
	return !rama_prepro_map_file_name_.empty();
}

/// @brief Set the names of the mainchain torsion potential maps to use to "".
///
void
ResidueTypeBase::reset_mainchain_torsion_potential_names() {
	rama_prepro_mainchain_torsion_potential_name_.clear();
	rama_prepro_map_file_name_.clear();
	rama_prepro_mainchain_torsion_potential_name_beforeproline_.clear();
	rama_prepro_map_file_name_beforeproline_.clear();
}

////////////////////////////////////////////////////////////////////////////
// adduct methods

void
ResidueTypeBase::report_adducts() const
{
	if ( defined_adducts_.size() == 0 ) return;

	for ( Size ii = 1 ; ii <= defined_adducts_.size() ; ++ii ) {
		Adduct const & add( defined_adducts_[ii] );
		TR.Debug << "Residue: " << name3() << " Adduct: " << add.adduct_name() <<
			" Atom name: " << add.atom_name() << std::endl;
	}
}

} // chemical
} // core


#ifdef    SERIALIZATION

template< class Archive >
void
core::chemical::ResidueTypeBase::save( Archive & arc ) const {
	using namespace core::chemical;

	arc( CEREAL_NVP( mode_ ) ); // enum core::chemical::TypeSetMode
	arc( CEREAL_NVP( atom_types_ ) ); // AtomTypeSetCOP
	arc( CEREAL_NVP( elements_ ) ); // ElementSetCOP
	arc( CEREAL_NVP( mm_atom_types_ ) ); // MMAtomTypeSetCOP
	arc( CEREAL_NVP( gasteiger_atom_types_ ) ); // gasteiger::GasteigerAtomTypeSetCOP
	arc( CEREAL_NVP( orbital_types_ ) ); // orbitals::OrbitalTypeSetCOP

	arc( CEREAL_NVP( aa_ ) ); // enum core::chemical::AA
	arc( CEREAL_NVP( backbone_aa_ ) ); // enum core::chemical::AA
	arc( CEREAL_NVP( na_analogue_ ) ); // enum core::chemical::AA
	arc( CEREAL_NVP( base_analogue_ ) ); // enum core::chemical::AA
	arc( CEREAL_NVP( base_name_ ) ); // std::string
	arc( CEREAL_NVP( name_ ) ); // std::string
	arc( CEREAL_NVP( name3_ ) ); // std::string
	arc( CEREAL_NVP( name1_ ) ); // char
	arc( CEREAL_NVP( interchangeability_group_ ) ); // std::string

	arc( CEREAL_NVP( rama_prepro_mainchain_torsion_potential_name_ ) ); // std::string
	arc( CEREAL_NVP( rama_prepro_mainchain_torsion_potential_name_beforeproline_ ) ); // std::string
	arc( CEREAL_NVP( rama_prepro_map_file_name_ ) );// std::string
	arc( CEREAL_NVP( rama_prepro_map_file_name_beforeproline_ ) );// std::string

	arc( CEREAL_NVP( net_formal_charge_ ) ); //signed long int
	arc( CEREAL_NVP( force_nbr_atom_orient_ ) ); // _Bool
	arc( CEREAL_NVP( remap_pdb_atom_names_ ) ); // _Bool
	arc( CEREAL_NVP( is_metapatched_ ) ); //bool
	arc( CEREAL_NVP( properties_ ) ); // ResiduePropertiesOP
	arc( CEREAL_NVP( rotamer_library_specification_ ) ); // rotamers::RotamerLibrarySpecificationOP

	arc( CEREAL_NVP( metal_binding_atoms_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( disulfide_atom_name_ ) ); // std::string
	arc( CEREAL_NVP( atom_aliases_ ) ); // std::map<std::string, std::string>
	arc( CEREAL_NVP( canonical_atom_aliases_ ) ); // std::map<std::string, std::string>
	arc( CEREAL_NVP( defined_adducts_ ) ); // utility::vector1<Adduct>

	arc( CEREAL_NVP( orbitals_ ) ); // utility::vector1<Orbital>
	arc( CEREAL_NVP( orbitals_index_ ) ); // std::map<std::string, int>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ResidueTypeBase::load( Archive & arc ) {
	using namespace core::chemical;

	arc( mode_ ); // enum core::chemical::TypeSetMode
	arc( atom_types_ ); // AtomTypeSetCOP
	arc( elements_ ); // ElementSetCOP
	arc( mm_atom_types_ ); // MMAtomTypeSetCOP
	arc( gasteiger_atom_types_ ); // gasteiger::GasteigerAtomTypeSetCOP
	arc( orbital_types_ ); // orbitals::OrbitalTypeSetCOP

	arc( aa_ ); // enum core::chemical::AA
	arc( backbone_aa_ ); // enum core::chemical::AA
	arc( na_analogue_ ); // enum core::chemical::AA
	arc( base_analogue_ ); // enum core::chemical::AA
	arc( base_name_ ); // std::string
	arc( name_ ); // std::string
	arc( name3_ ); // std::string
	arc( name1_ ); // char
	arc( interchangeability_group_ ); // std::string

	arc( rama_prepro_mainchain_torsion_potential_name_ ); // std::string
	arc( rama_prepro_mainchain_torsion_potential_name_beforeproline_ ); // std::string
	arc( rama_prepro_map_file_name_ ); // std::string
	arc( rama_prepro_map_file_name_beforeproline_ ); // std::string

	arc( net_formal_charge_ ); //signed long int
	arc( force_nbr_atom_orient_ ); // _Bool
	arc( remap_pdb_atom_names_ ); // _Bool
	arc( is_metapatched_ ); //bool
	arc( properties_ );
	arc( rotamer_library_specification_ );

	arc( metal_binding_atoms_ ); // utility::vector1<std::string>
	arc( disulfide_atom_name_ ); // std::string
	arc( atom_aliases_ ); // std::map<std::string, std::string>
	arc( canonical_atom_aliases_ ); // std::map<std::string, std::string>
	arc( defined_adducts_ ); // utility::vector1<Adduct>

	arc( orbitals_ ); // utility::vector1<Orbital>
	arc( orbitals_index_ ); // std::map<std::string, int>
}


SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ResidueTypeBase );
CEREAL_REGISTER_TYPE( core::chemical::ResidueTypeBase )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_ResidueTypeBase )
#endif // SERIALIZATION
