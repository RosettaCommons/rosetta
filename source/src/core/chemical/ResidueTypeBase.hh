// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ResidueTypeBase.hh
/// @brief Method declarations and simple accessors/getters for ResidueTypeBase
/// @author
/// Phil Bradley
/// Rocco Moretti (rmorettiase@gmail.com)
/// Steven Combs
/// Vikram K. Mulligan - properties for D-, beta- and other noncanonicals
/// Jason W. Labonte (code related to properties, lipids, carbohydrates, and other non-AAs)

#ifndef INCLUDED_core_chemical_ResidueTypeBase_hh
#define INCLUDED_core_chemical_ResidueTypeBase_hh

// Unit headers
#include <core/chemical/ResidueTypeBase.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh> // For ResidueTypeCOP in base type return

// Package headers
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.fwd.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/Orbital.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>
#include <utility/pointer/deep_copy.hh>

// C++ headers
#include <map>

#ifdef    SERIALIZATION
#include <utility/serialization/serialization.hh>
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

/// @brief A base class for definiting types of residues
/// @details
/// This class is a base class for the "chemical" information for residues in Rosetta.
/// It's a base class for two different classes. ResidueType, which is the main class used by Rosetta,
/// but which is a fixed-content type, and MutableResidueType, which is a modifiable class.
/// The two type also differ also in how Atoms information is represented.
/// The atom information in MutableResidueType is represented in a molecular graph, which is more convenient for modification,
/// with atoms referenced primarily by "vertex descriptor" in the graph.
/// The atom information in ResidueType is in a struct-of-arrays format, which allows better optimization, with atoms
/// in a defined order.
/// There are also differences in guarantees for each class. ResidueType has certain atom ordering guarantees which MutableResidueType lacks.
///
/// To accomodate the two different classes, this base class only stores the common, (non-atom dependent) information,
/// along with some convenient accessors.
///
/// See the documentation of each class for more information.
///
/// A note on things which belong in this class: Only include things which are atom-independent,
/// and are "primary" information. (Things which can/should be derived from atoms or other info
/// should be placed in the plain ResidueType class,
/// and updated on MutableResidueType -> ResidueType conversion.
/// You can also include (as a virtual method) functions needed for PatchSelector function,
/// but any atom indexing/references should be name-based.
class ResidueTypeBase : public utility::VirtualBase
{
protected:
	ResidueTypeBase(); // protected, not deleted, as serialization needs it.

	/// @brief constructor
	/// @details We use the AtomTypeSet object to assign atom_types to atoms inside add_atom,
	/// and to identify (polar) hydrogens, acceptors, etc.
	/// This is protected as we don't want to instantiate a plain ResidueType directly.
	ResidueTypeBase(
		AtomTypeSetCOP atom_types,
		ElementSetCOP element_types,
		MMAtomTypeSetCOP mm_atom_types,
		orbitals::OrbitalTypeSetCOP orbital_types
	);

protected: // Don't allow general copying of a ResidueTypeBase

	ResidueTypeBase( ResidueTypeBase const & ) = default;

	/// @brief Copies  <src>  into the ResidueTypeBase
	ResidueTypeBase &
	operator=( ResidueTypeBase const & ) = default;

public:

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Type Set Funcitons          ////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	TypeSetMode
	mode() const { return mode_; }

	/// @brief access by reference the atomset for which this residue is constructed
	AtomTypeSet const &
	atom_type_set() const
	{
		debug_assert( atom_types_ );
		return *atom_types_;
	}

	/// @brief access by const pointer the atomset for which this residue is constructed
	AtomTypeSetCOP
	atom_type_set_ptr() const
	{
		return atom_types_;
	}

protected: // deep magic functions - don't expose publically

	void atom_type_set( AtomTypeSetCOP setting );

public:

	/// @brief access by reference the element set for which this residue is constructed
	ElementSet const &
	element_set() const
	{
		debug_assert( elements_ );
		return *elements_;
	}

	/// @brief access by const pointer the element set for which this residue is constructed
	ElementSetCOP
	element_set_ptr() const
	{
		return elements_;
	}

	MMAtomTypeSetCOP
	mm_atom_types_ptr() const {
		return mm_atom_types_;
	}

	gasteiger::GasteigerAtomTypeSetCOP
	gasteiger_atom_typeset() const;

	void
	set_gasteiger_atom_typeset( gasteiger::GasteigerAtomTypeSetCOP setting );

	/// @brief Get the MM atom_type for this atom by its index number in this residue
	orbitals::OrbitalTypeSetCOP
	orbital_types_ptr() const {
		return orbital_types_;
	}

	void
	set_orbital_typeset( orbitals::OrbitalTypeSetCOP setting ) { orbital_types_ = setting; }

	/// @brief Does this residue type define orbital types?
	bool has_orbital_types() const { return orbital_types_ != nullptr; }

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Summary Functions              /////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief number of atoms
	virtual
	Size
	natoms() const = 0;

	/// @brief number of heavy atoms
	virtual
	Size
	nheavyatoms() const = 0;

	/// @brief number of bonds
	virtual
	Size
	nbonds() const = 0;

	/// @brief Counts the number of virtual atoms and returns the count.
	/// @details The virtual count is not stored in the resiude type.  This count is performed on the fly, and
	///can hurt performance if reapeatedly carried out.  Not intended for use in large loops -- instead, call
	///once and store the value.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	virtual
	Size
	n_virtual_atoms () const = 0;

	/// @brief is this atom present in this residue?
	virtual
	bool
	has( std::string const & atom_name ) const = 0;

#ifdef WIN32  // Fixes incorrect cast on WIN32 where has("string") actually calls has( VD )
	inline bool has( const char *name ) const
	{
		return has( std::string(name) );
	}
#endif

	/// @brief Get the nominal net formal charge on this residue type.
	/// @details This may not match the sum of the formal charges on the atoms
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	signed long int net_formal_charge() const { return net_formal_charge_; }

	/// @brief Set the nominal net formal charge on this residue type.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void net_formal_charge( signed long int charge_in ) {
		net_formal_charge_ = charge_in;
	}

	virtual
	void
	show_all_atom_names( std::ostream & out ) const = 0;

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Data Accessors                 /////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	utility::vector1< std::string > const &
	get_metal_binding_atoms() const { return metal_binding_atoms_; }

	/// @brief Gets disulfide atom name
	/// @author Andrew M. Watkins (amw579@nyu.edu).
	std::string const &
	get_disulfide_atom_name() const {
		return disulfide_atom_name_;
	}

	/// @brief Sets disulfide atom name
	/// @author Andrew M. Watkins (amw579@nyu.edu).
	void
	set_disulfide_atom_name( std::string const & n ) {
		disulfide_atom_name_ = n;
	}

	/// @brief returns atom aliases
	std::map< std::string, std::string > const &
	atom_aliases() const {
		return atom_aliases_;
	}

	/// @brief returns atom alias
	std::string const &
	atom_alias( std::string const & name) const {
		debug_assert( atom_aliases_.count(name) == 1 );
		return atom_aliases_.at(name);
	}

	/// @brief returns atom aliases
	std::map< std::string, std::string > const &
	canonical_atom_aliases() const {
		return canonical_atom_aliases_;
	}

	/// @brief returns atom alias
	std::string const &
	canonical_atom_alias( std::string const & name) const {
		debug_assert( canonical_atom_aliases_.count(name) == 1 );
		return canonical_atom_aliases_.at(name);
	}

protected: // Non-const accessors

	utility::vector1< std::string > &
	metal_binding_atoms() { return metal_binding_atoms_; }

	/// @brief returns atom aliases
	std::map< std::string, std::string > &
	atom_aliases() {
		return atom_aliases_;
	}

	/// @brief returns atom aliases
	std::map< std::string, std::string > &
	canonical_atom_aliases() {
		return canonical_atom_aliases_;
	}

public:

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	////////////////          Orbital Functions     //////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	Orbital const & orbital(Size const orbital_index) const;
	Orbital const & orbital(std::string const & orbital_name) const;

	orbitals::OrbitalType const &
	orbital_type(Size const orbital_index) const;

	/// @brief number of orbitals
	Size
	n_orbitals() const;

	/// @brief is this orbital present in this residue?
	bool
	has_orbital( std::string const & orbital_name ) const
	{
		return ( orbitals_index_.find( orbital_name ) != orbitals_index_.end() );
	}

	/// @brief get orbital index by name
	core::Size
	orbital_index( std::string const & name ) const;

protected: // Modifiable accessors

	utility::vector1< Orbital > &
	get_orbitals() { return orbitals_; }

	std::map< std::string, int > &
	get_orbitals_index() { return orbitals_index_; }

public:

	//////////////////////////////////////////////////////////////////////
	///////////////// properties /////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief Access the collection of properties for this ResidueTypeBase.
	ResidueProperties const & properties() const { return *properties_; };

	/// @brief Set the collection of properties for this ResidueTypeBase.
	void set_properties( ResiduePropertiesOP properties );

	/// @brief Add a property to this ResidueTypeBase.
	void add_property( std::string const & property );

	/// @brief Add a property to this ResidueType, by properties enum.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void add_property( core::chemical::ResidueProperty const property );

	void set_adduct_flag( bool adduct_in );

	/// @brief Add a numeric property.
	void add_numeric_property( std::string const & tag, core::Real value );

	/// @brief Add a string property.
	void add_string_property( std::string const & tag, std::string value );

	/// @brief Delete a property of this ResidueType.
	void delete_property( std::string const & property );

	/// @brief Delete a property of this ResidueType, by properties enum.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void delete_property( core::chemical::ResidueProperty const property );

	/// @brief Is this ResidueTypeBase a base type?
	virtual
	bool is_base_type() const = 0;

	/// @brief Get a pointer to this ResidueTypeBase's base ResidueTypeBase.
	/// @details NOTE: Behavior when `is_base_type() == true` varies by subclass!
	virtual
	ResidueTypeCOP get_base_type_cop() const = 0;

	/// @brief Set that this is a metapatched ResidueType.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void set_metapatched() { is_metapatched_ = true; }

	/// @brief Get whether this is a metapatched ResidueType.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	inline bool is_metapatched() const { return is_metapatched_; }

	// Note: most of the is_* functions are deliberately only on plain ResidueType.
	// These depend on the property settings, which may not be consistent during MutableResidueType building.
	// The following are only a small subset, useful during MutableResidueType construction itself.

	/// @brief is polymer?
	bool is_polymer() const { return properties().has_property( POLYMER ); }

	/// @brief is protein?
	bool is_protein() const { return properties().has_property( PROTEIN ); }

	/// @brief is DNA?
	bool is_DNA() const { return properties().has_property( DNA ); }

	/// @brief is RNA?
	bool is_RNA() const { return properties().has_property( RNA ); }

	/// @brief is this a d-amino acid?
	bool is_d_aa() const { return properties().has_property( D_AA ); }

	/// @brief is this an l-amino acid?
	bool is_l_aa() const { return properties().has_property( L_AA ); }

	/// @brief is this a d-RNA?
	bool is_d_rna() const { return properties().has_property( D_RNA ); }

	/// @brief is this an l-RNA?
	bool is_l_rna() const { return properties().has_property( L_RNA ); }

	/// @brief is peptoid?
	bool is_peptoid() const { return properties().has_property( PEPTOID ); }

	/// @brief Is this a peptoid with a chiral side-chain with "R" chirality?
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool is_r_peptoid() const { return properties().has_property( R_PEPTOID ); }

	/// @brief Is this a peptoid with a chiral side-chain with "S" chirality?
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool is_s_peptoid() const { return properties().has_property( S_PEPTOID ); }

	/// @brief is this an achiral backbone?
	bool is_achiral_backbone() const { return properties().has_property( ACHIRAL_BACKBONE ); }

	/// @brief Does this have an achiral sidechain?
	/// @details Includes gly and aib, and most (but not all) peptoids.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool is_achiral_sidechain() const { return properties().has_property( ACHIRAL_SIDECHAIN ); }

	/// @brief  Generic property access.
	bool has_property( std::string const & property ) const;

	/// @brief  Generic property access, by ResidueProperty.
	bool has_property( ResidueProperty const property ) const;

	/// @brief Get a numeric property, if it exists.
	core::Real get_numeric_property(std::string const & tag) const;

	/// @brief Get a string property, if it exists.
	std::string get_string_property(std::string const & tag) const;

	/// @brief Add a variant type to this ResidueTypeBase.
	void add_variant_type( VariantType const variant_type );

	/// @brief Add a variant type to this ResidueTypeBase by string.
	void add_variant_type( std::string const & variant_type );

	/// @brief Remove a variant type to this ResidueTypeBase.
	void remove_variant_type( VariantType const variant_type );

	/// @brief Remove a variant type to this ResidueTypeBase by string.
	void remove_variant_type( std::string const & variant_type );

	/// @brief  Generic variant access.
	bool has_variant_type( VariantType const variant_type ) const;

	// TODO: Find a way to remove this; it only exists because of how ResidueTypeSelectors are currently written. ~Labonte
	/// @brief  Generic variant access by string.
	bool has_variant_type( std::string const & variant_type ) const;

	/// @brief  Turn on the ability to create VariantTypes "on-the-fly".
	void enable_custom_variant_types();

	/// @brief get all the variant types for this ResidueTypeBase
	/// @details This will include both on-the-fly custom variants defined by string AND string equivalents
	/// of standard, enumerated variants.
	///         -- rhiju (merging roccomoretti/restypeset_fiddle)
	utility::vector1< std::string >
	variant_types() const;

	/// @brief Get a vector of VariantType enums for this ResidueTypeBase.
	/// @details This ONLY includes standard, enum-based variants, not on-the-fly custom variants.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	utility::vector1< VariantType > variant_type_enums() const;

	/// @brief Get a list of custom VariantType strings for this ResidueTypeBase (by const reference).
	/// @details This ONLY includes custom, on-the-fly variants, not standard variants.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	utility::vector1< std::string > const & custom_variant_types() const;


	/// @brief set our aa-type (could be "UNK")
	void
	aa( AA const & type )
	{
		aa_ = type;
	}
	/// @brief set our aa-type (could be "UNK")
	void
	aa( std::string const & type )
	{
		aa_ = aa_from_name( type );
	}

	/// @brief AA to use for backbone scoring
	void
	backbone_aa( std::string const & type )
	{
		backbone_aa_ = aa_from_name( type );
	}

	/// @brief Set AA to use for backbone scoring directly (without string conversion).
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void
	backbone_aa( core::chemical::AA const type ) {
		backbone_aa_ = type;
	}

	/// @brief NA to use for base-specific generalization (can be more
	/// forgiving than na_analogue for new NA backbones)
	void
	base_analogue( std::string const & type )
	{
		base_analogue_ = aa_from_name( type );
	}

	/// @brief NA to use for fragment sampling and some scoring purposes
	void
	na_analogue( std::string const & type )
	{
		na_analogue_ = aa_from_name( type );
	}

	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	///////////////// Names     /////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////

	/// @brief Get this ResidueTypeBase's base name (shared with other residue types with the same base type).
	///
	std::string const &
	base_name() const {
		return base_name_;
	}

	/// @brief get our (unique) residue name
	std::string const &
	name() const
	{
		return name_;
	}

	/// @brief Set this ResidueTypeBase's base name (shared with other residue types with the same base type).
	///
	void base_name( std::string const &base_name_in ) {
		base_name_ = base_name_in;
	}

	/// @brief set our (unique) residue name
	void
	name( std::string const & name_in )
	{
		name_ = name_in;
		properties_->set_parent_name( name_in );
	}


	/// @brief get our 3letter code.  This is set in the
	/// ResidueTypeBase .params file through the IO_STRING
	/// tag along with the name1 string
	/// NOTE: The "name3" is not necessarily three characters
	/// long. e.g. Metal ions may be only two characters.
	/// @details If you need three characters, the PDB convention
	/// is to right pad.
	std::string const &
	name3() const
	{
		return name3_;
	}

	/// @brief set our 3letter code
	void
	name3( std::string const & name_in )
	{
		name3_ = name_in;
	}


	/// @brief get our 1letter code.  This is set in the
	/// ResidueTypeBase .params file through the IO_STRING
	/// tag along with the name3 string.
	char
	name1() const
	{
		return name1_;
	}


	/// @brief set our 1letter code
	void
	name1( char const code )
	{
		name1_ = code;
	}

	/// @brief get our interchangeability-group id.  Used to
	/// determine if two residue types are equivalent, except
	/// for their variant status.  E.g. ResidueTypeBases ALA and
	/// ALA_Nterm would be part of the same interchangeability
	/// group.  This has a degree of subjectivity; are TYR and
	/// pTYR in the same interchangeability group?  Probably not.
	/// This data can be specified in the ResidueTypeBases .params
	/// file with the INTERCHANGEABILITY_GROUP tag.
	std::string
	interchangeability_group() const {
		return interchangeability_group_;
	}

	/// @brief set our interchangeability-group id
	void
	interchangeability_group( std::string const & setting ) {
		interchangeability_group_ = setting;
	}

	/// @brief Turn on geometry-based atom renaming when loading this residue type from PDB files
	void
	remap_pdb_atom_names( bool rename )
	{
		remap_pdb_atom_names_ = rename;
	}

	/// @brief Are we using geometry-based atom renaming when loading this residue type from PDB
	bool remap_pdb_atom_names() const
	{
		return remap_pdb_atom_names_;
	}

	/// @brief our traditional residue type, if any
	///
	/// @details Used for knowledge-based scores, dunbrack, etc.
	/// could be "aa_unk".
	///
	/// AA is an enum.  There are values for the 20 standard amino acids, the
	/// 19 canonical D-amino acids, common beta-amino acids and nucleic acids,
	/// and aa_unk as a general catch-all.
	AA const &
	aa() const
	{
		return aa_;
	}

	/// @brief Returns the amino acid type to be used for backbone scoring (rama and p_aa_pp).
	AA const &
	backbone_aa() const
	{
		if ( backbone_aa_==aa_unk ) return aa_;
		return backbone_aa_;
	}

	/// @brief Returns the nucleic acid type to be used for base features
	AA const &
	base_analogue() const
	{
		if ( base_analogue_ == aa_unp ) return aa_;
		return base_analogue_;
	}

	/// @brief Returns the nucleic acid type to be used for fragment sampling/scoring.
	AA const &
	na_analogue() const
	{
		if ( na_analogue_ == aa_unp ) return aa_;
		return na_analogue_;
	}

protected: // The raw versions should only be used internally by derived classes.

	AA const &
	backbone_aa_raw() const { return backbone_aa_; }

	AA const &
	base_analogue_raw() const { return base_analogue_; }

	AA const &
	na_analogue_raw() const { return na_analogue_; }

public:

	////////////////////////////////////////////////////////////
	///////////// Rotamer library //////////////////////////////
	////////////////////////////////////////////////////////////

	void
	rotamer_library_specification( rotamers::RotamerLibrarySpecificationOP );

	rotamers::RotamerLibrarySpecificationCOP
	rotamer_library_specification() const;

	/// @brief Nonconst access to the RotamerLibrarySpecification.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	rotamers::RotamerLibrarySpecificationOP
	rotamer_library_specification_nonconst();

	/// @brief Remove any rotamer library specifications attached to this ResidueTypeBase.
	/// @details After this operation, the rotamer_library_specification() method returns a NULL pointer.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void strip_rotamer_library_specification();

	////////////////////////////////////////////////////////////////////////////
	/// dihedral methods
public:

	/// @brief Return force_nbr_atom_orient_, used to control orient atoms selected by select_orient_atoms
	bool force_nbr_atom_orient() const
	{
		return force_nbr_atom_orient_;
	}

	/// @brief Set force_nbr_atom_orient_, used to control orient atoms selected by select_orient_atoms
	void
	force_nbr_atom_orient( bool force_orient )
	{
		force_nbr_atom_orient_ = force_orient;
	}

	////////////////////////////////////////////////////////////////////////////
	// rama methods
	////////////////////////////////////////////////////////////////////////////

	/// @brief Get the key name for the mainchain torsion potential used by the RamaPrePro score term.
	/// @details Stored internally as a string for base residue types.  Empty string is stored by default for derived
	/// residue types (in which case this function returns the string stored in the base ResidueTypeBase), though this can be overridden.
	/// @note Different maps are used for preproline positions and non-preproline positions.  The boolean determines which map
	/// we're interested in.
	std::string const & get_rama_prepro_mainchain_torsion_potential_name( bool const pre_proline_position ) const;

	/// @brief Do the rama_prepro mainchain torsion potentials of this residue match another?
	///
	bool mainchain_potentials_match( ResidueTypeBase const &other ) const;

	/// @brief Set the key name for the mainchain torsion potential used by the RamaPrePro score term.
	/// @details Stored internally as a string for base residue types.  Empty string is stored by default for derived
	/// residue types (pointing the function to the base type), though this can be overridden using this function.
	/// @note Different maps are used for preproline positions and non-preproline positions.  The boolean determines which map
	/// we're interested in.
	void set_rama_prepro_mainchain_torsion_potential_name( std::string const &name_in, bool const pre_proline_position);

	/// @brief Get the file name for the mainchain torsion potential used by the RamaPrePro score term.
	/// @details Stored internally as a string for base residue types.  Empty string is stored by default for derived
	/// residue types (in which case this function returns the string stored in the base ResidueTypeBase), though this can be overridden.
	/// @note Different maps are used for preproline positions and non-preproline positions.  The boolean determines which map
	/// we're interested in.
	std::string const & get_rama_prepro_map_file_name(bool const pre_proline_position) const;

	/// @brief Set the file name for the mainchain torsion potential used by the RamaPrePro score term.
	/// @details Stored internally as a string for base residue types.  Empty string is stored by default for derived
	/// residue types (pointing the function to the base type), though this can be overridden using this function.
	/// @note Different maps are used for preproline positions and non-preproline positions.  The boolean determines which map
	/// we're interested in.
	void set_rama_prepro_map_file_name( std::string const &filename_in, bool const pre_proline_position );

	/// @brief Returns true if and only if (a) this is not a base type, AND (b) there is a rama_prepro_mainchain_torsion_map_file_name_
	/// defined for this ResidueTypeBase (which is presumably different from that of the base type).
	/// @details If pre_proline_position is true, checks rama_prepro_mainchain_torsion_map_file_name_beforeproline_ instead of
	/// rama_prepro_mainchain_torsion_potential_name_.
	bool defines_custom_rama_prepro_map( bool const pre_proline_position ) const;

	/// @brief Set the names of the mainchain torsion potential maps to use to "".
	/// @details Also resets the mainchain torsion potential filename strings.
	void reset_mainchain_torsion_potential_names();

	/// @brief  Generate string representation of ResidueTypeBase for debugging purposes.
	virtual
	void show( std::ostream & output=std::cout, bool output_atomic_details=false ) const = 0;

public:
	////////////////////////////////////////////////////////////////////////////
	// adduct methods
	/// @brief get the adducts defined for this residue
	utility::vector1< Adduct > const &
	defined_adducts() const
	{
		return defined_adducts_;
	}

	void
	add_adduct( Adduct & adduct_in )
	{
		defined_adducts_.push_back( adduct_in );
	}

	void
	report_adducts() const;

private:

	//////////////////////////////////////////////////////////////////
	// General Typeset information

	/// @brief The TypeSet mode for this ResidueTypeBase
	TypeSetMode mode_ = INVALID_t;

	/// @brief The type set for Rosetta Atom types.
	/// Required to be non-null for standard ResidueTypes
	AtomTypeSetCOP atom_types_;

	/// @brief The set for element objects. -- can be null.
	ElementSetCOP elements_;

	/// @brief The set for MMAtomTypes. -- can be null.
	MMAtomTypeSetCOP mm_atom_types_;

	/// @brief The set for GasteigerAtomTypes. -- can be null.
	gasteiger::GasteigerAtomTypeSetCOP gasteiger_atom_types_;

	/// @brief The set for OrbitalTypes. -- can be null.
	orbitals::OrbitalTypeSetCOP orbital_types_;

	///////////////////////////////////////////////////////////////////////////
	// General Residue properties

	/// @brief standard rosetta aa-type for knowledge-based potentials, may be aa_unk
	// aa_ = THIS residue's aa-type
	// rotamer_aa_ = the aa-type on which rotamers will be based
	// backbone_aa_ = the aa-type on which the backbone scoring (rama, p_aa_pp) will be based
	// na_analogue_ = the aa-type for a nucleic acid (generally, canonical) that
	// fragment sampling will use.
	// base_analogue_ = the aa-type for a nucleic acid for base-specific features
	// like WC base pair atoms etc.
	AA aa_ = aa_unk, backbone_aa_ = aa_unk, na_analogue_ = aa_unp, base_analogue_ = aa_unp;

	/// @brief Residue id for the base type (i.e. sans variants).
	/// @details Does not accumulate VariantType names (e.g. "ALA").
	std::string base_name_;

	/// @brief Residue id -- Should be unique within a ResidueTypeSet
	/// @details Accumulates VariantType names (e.g. "ALA:NtermProteinFull").
	std::string name_;

	/// @brief PDB-file id, need not be unique
	/// NOTE: The "name3" is not necessarily three characters long.
	/// e.g. Metal ions may be only two characters.
	std::string name3_;

	/// @brief one-letter code, also not necessarily unique
	char name1_ = ' ';

	/// @brief interchangeability group lets a ResidueTypeBase claim to be functionally
	/// interchangeable with any other ResidueTypeBase in the same group.  This
	/// is used by the packer to decide which ResidueTypeBase from a desired group
	/// has the right set of variants to be placed at a particular position.
	/// E.g. if the interchangeability group is "ALA" and the packer is building
	/// rotamers for residue 1, (the N-terminal residue) then, the packer will
	/// select the "ALA:NTermProteinFull" ResidueTypeBase and build rotamers for it.
	std::string interchangeability_group_;

	/// @brief The name of the mainchain torsion potential used by this ResidueTypeBase for the RamaPrePro score term.
	/// @details If blank, the base type's mainchain torsion potential is used.
	/// @note This one is for maps used if this residue is not immediately before a proline.
	std::string rama_prepro_mainchain_torsion_potential_name_;

	/// @brief The name of the mainchain torsion potential used by this ResidueTypeBase for the RamaPrePro score term.
	/// @details If blank, the base type's mainchain torsion potential is used.
	/// @note This one is for maps used if this residue occurs before a proline.
	std::string rama_prepro_mainchain_torsion_potential_name_beforeproline_;

	/// @brief The filename of the mainchain torsion potential used by this ResidueTypeBase for the RamaPrePro score term.
	/// @details If blank, the base type's file is used.  If that's blank too, then this ResidueTypeBase isn't scored by RamaPrePro.
	/// @note This one is for maps used if this residue is not immediately before a proline.
	std::string rama_prepro_map_file_name_;

	/// @brief The filename of the mainchain torsion potential used by this ResidueTypeBase for the RamaPrePro score term.
	/// @details If blank, the base type's file is used.  If that's blank too, then this ResidueTypeBase isn't scored by RamaPrePro.
	/// @note This one is for maps used if this residue occurs before a proline.
	std::string rama_prepro_map_file_name_beforeproline_;

	/// @brief What is the total formal charge of the ResidueType?
	/// (Note that this may be different from the sum of the formal charges on each atom.)
	signed long int net_formal_charge_ = 0;

	/// Controls which atoms are selected by "select_orient_atoms",
	/// used to overlay residues during packing.
	bool force_nbr_atom_orient_ = false;

	/// @brief Should we attempt to rename atoms for this residue type
	/// when we read in PDB files?
	bool remap_pdb_atom_names_ = false;

	/// @brief Residue properties as defined in the residue topology (.params) files
	utility::pointer::DeepCopyOP< ResidueProperties > properties_;

	/// @brief Is this ResidueType the product of metapatching?
	/// @details False by default; true if and only if a metapatch contributed to this ResidueType.
	bool is_metapatched_ = false;

	/// @brief How to construct a rotamer library for this ResidueTypeBase
	utility::pointer::DeepCopyOP< rotamers::RotamerLibrarySpecification > rotamer_library_specification_;

	///////////////////////////////////////////////////////////////////////////
	// Atom-name-based atom annotations

	/// @brief Names of all of the atoms that are able to make a bond to a metal, for metal-binding residue types
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	utility::vector1 < std::string > metal_binding_atoms_;

	/// @brief Name of the disulfide-forming atom, if any
	std::string disulfide_atom_name_ = "NONE";

	/// @brief A mapping of alias atom names to canonical atom names
	/// Will be added to atom_name_to_vd_ during finalization
	std::map< std::string, std::string > atom_aliases_;

	/// @brief A map of canonical atom names to atom aliases (white space included)
	std::map< std::string, std::string > canonical_atom_aliases_;

	/// @brief Adducts defined for this residue
	utility::vector1< Adduct > defined_adducts_;

	/// @brief The orbitals on the ResidueTypeBase, if any.
	utility::vector1< Orbital > orbitals_;

	/// @brief index lookup for orbitals based on atom name
	/// Updated in add_orbital()
	std::map< std::string, int > orbitals_index_;

public:

#ifdef    SERIALIZATION
public:
	friend class cereal::access;
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};  // class ResidueTypeBase

} // chemical
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_ResidueTypeBase )
#endif // SERIALIZATION

#endif // INCLUDED_core_chemical_Residues_HH
