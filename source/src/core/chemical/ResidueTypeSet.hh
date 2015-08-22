// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/ResidueTypeSet.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_ResidueTypeSet_hh
#define INCLUDED_core_chemical_ResidueTypeSet_hh


// Unit headers
#include <core/chemical/ResidueTypeSet.fwd.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSelector.fwd.hh>

// STL headers
#include <list>

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Patch.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>
#include <set>

//Auto Headers
#include <core/chemical/Adduct.fwd.hh>

#ifdef WIN32
#include <core/chemical/Adduct.hh>
#endif

namespace core {
namespace chemical {

///  A collection of ResidueType defined
/**
One thing that is not nailed down is whether a single ResidueSet can have ResidueType's with
different AtomTypeSets. I've left open this possibility currently although there isnt any code
for this yet (PB-07/07)
**/


class ResidueTypeSet : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< ResidueTypeSet >
{
public:
	typedef std::list< AA >::const_iterator AAsIter;
	typedef std::map< std::string, ResidueTypeCOP >::const_iterator const_residue_iterator;

public:


	/// @brief default c-tor
	ResidueTypeSet();

	/// @brief constructor from directory
	ResidueTypeSet(
		std::string const & name,
		std::string const & directory
	);

	virtual ~ResidueTypeSet();

	void init(
		std::vector< std::string > const & extra_res_param_files = std::vector< std::string >(),
		std::vector< std::string > const & extra_patch_files = std::vector< std::string >()
	);

	/// self pointers
	inline ResidueTypeSetCOP get_self_ptr() const { return shared_from_this(); }
	inline ResidueTypeSetOP  get_self_ptr() { return shared_from_this(); }
	inline ResidueTypeSetCAP get_self_weak_ptr() const { return ResidueTypeSetCAP( shared_from_this() ); }
	//inline ResidueTypeSetAP  get_self_weak_ptr() { return ResidueTypeSetAP( shared_from_this() ); }

	/// @brief name of the residue type set
	std::string const &
	name() const {
		return name_;
	}

	AtomTypeSetCOP atom_type_set() const { return atom_types_; }
	ElementSetCOP element_set() const { return elements_; }
	MMAtomTypeSetCOP mm_atom_type_set() const { return mm_atom_types_; }
	orbitals::OrbitalTypeSetCOP orbital_type_set() const { return orbital_types_; }

	void atom_type_set(AtomTypeSetCOP atom_types) {
		runtime_assert( !atom_types_ ); // Don't change a set default.
		atom_types_ = atom_types;
	}
	void element_set(ElementSetCOP elements) {
		runtime_assert( !elements_ ); // Don't change a set default.
		elements_ = elements;
	}
	void mm_atom_type_set(MMAtomTypeSetCOP mm_atom_types) {
		runtime_assert( !mm_atom_types_ ); // Don't change a set default.
		mm_atom_types_ = mm_atom_types;
	}
	void orbital_type_set(orbitals::OrbitalTypeSetCOP orbital_types) {
		runtime_assert( !orbital_types_ ); // Don't change a set default.
		orbital_types_ = orbital_types;
	}

	/// @brief read a list of residue types
	void
	read_list_of_residues(
		std::string const & list_filename
	);

	void
	read_files(
		utility::vector1< std::string > const & filenames
	);

	void
	apply_patches(
		utility::vector1< std::string > const & filenames
	);

	/// @brief apply patches to base ResidueType to generate variant ResidueTyes
	void
	place_adducts();

	/// @brief adds a new residue type to the set
	void
	add_residue_type( ResidueTypeOP new_type );

	void
	add_residue_type( std::string const &  filename );

	/// @brief delete a residue type from the set (Use with care)
	void
	remove_residue_type(std::string const & name);

	/// @brief Create correct combinations of adducts for a residue type
	void create_adduct_combinations(
		ResidueType const & rsd,
		std::map< std::string, int > ref_map,
		std::map< std::string, int > count_map,
		utility::vector1< bool > add_mask,
		utility::vector1< Adduct >::const_iterator work_iter
	);

	/// @brief query ResidueTypes by their interchangeability-group name
	ResidueTypeCOPs const &
	interchangeability_group_map_DO_NOT_USE( std::string const & name ) const;

	/// @brief Does this ResidueTypeSet have ResidueTypes with the given interchangability group?
	bool
	has_interchangeability_group( std::string const & name ) const;

	/// @brief query ResidueTypes by their 3-letter name
	ResidueTypeCOPs const &
	name3_map_DO_NOT_USE( std::string const & name ) const;

	/// @brief query ResidueType by its unique residue id.
	///
	/// @details since residue id is unique, it only returns
	/// one residue type or exit without match.
	virtual
	ResidueType const &
	name_map( std::string const & name ) const;

	/// @brief Get the base ResidueType with the given aa type
	/// @details Returns 0 if one does not exist.
	ResidueTypeCOP
	get_representative_type_aa( AA aa ) const;

	/// @brief Get the base ResidueType with the given name1
	/// @details Returns 0 if one does not exist.
	ResidueTypeCOP
	get_representative_type_name1( char name1 ) const;

	/// @brief Get the base ResidueType with the given name3
	/// @details Returns 0 if one does not exist.
	ResidueTypeCOP
	get_representative_type_name3( std::string const &  name3 ) const;

	/// @brief Get the base ResidueType with the given aa type and variants
	/// @details Returns 0 if one does not exist.
	/// The returned type will have at least all the variants given, but may have more
	/// if a minimal variant type isn't availible.
	ResidueTypeCOP
	get_representative_type_aa( AA aa, utility::vector1< std::string > const & variants ) const;

	/// @brief Get the base ResidueType with the given name1 and variants
	/// @details Returns 0 if one does not exist.
	/// The returned type will have at least all the variants given, but may have more
	/// if a minimal variant type isn't availible.
	ResidueTypeCOP
	get_representative_type_name1( char name1, utility::vector1< std::string > const & variants ) const;

	/// @brief Get the base ResidueType with the given name3 and variants
	/// @details Returns 0 if one does not exist.
	/// The returned type will have at least all the variants given, but may have more
	/// if a minimal variant type isn't availible.
	ResidueTypeCOP
	get_representative_type_name3( std::string const &  name3, utility::vector1< std::string > const & variants ) const;

	/// @brief Gets all non-patched types with the given aa type
	ResidueTypeCOPs
	get_base_types_aa( AA aa ) const;

	/// @brief Get all non-patched ResidueTypes with the given name1
	ResidueTypeCOPs
	get_base_types_name1( char name1 ) const;

	/// @brief Get all non-patched ResidueTypes with the given name3
	ResidueTypeCOPs
	get_base_types_name3( std::string const &  name3 ) const;

	/// @brief Check if a base type (like "SER") generates any types with another name3 (like "SEP")
	bool
	generates_patched_residue_type_with_name3( std::string const & base_residue_name, std::string const & name3 ) const;

	/// @brief Gets all types with the given aa type and variants
	/// @details The number of variants must match exactly.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	ResidueTypeCOPs
	get_all_types_with_variants_aa( AA aa, utility::vector1< std::string > const & variants ) const;

	/// @brief Get all non-patched ResidueTypes with the given name1
	/// @details The number of variants must match exactly.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	ResidueTypeCOPs
	get_all_types_with_variants_name1( char name1, utility::vector1< std::string > const & variants ) const;

	/// @brief Get all non-patched ResidueTypes with the given name3
	/// @details The number of variants must match exactly.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	ResidueTypeCOPs
	get_all_types_with_variants_name3( std::string const &  name3, utility::vector1< std::string > const & variants ) const;

	// IF YOU COMMENT THIS BACK IN, TELL ME WHY BY EMAIL: aleaverfay@gmail.com
	///// @brief query ResidueType by its unique residue id returning a non-const reference to that ResidueType
	///// (Is this really a good idea?  Who is using this functionality?  ResidueTypes are supposed
	///// to be constant once deposited inside a ResidueTypeSet or lots of downstream assumptions
	///// are invalidated!)
	//ResidueType &
	//nonconst_name_map( std::string const & name );

	/// @brief query if a ResidueType of the unique residue id (name) is present.
	bool has_name( std::string const & name ) const;

	/// @brief query if any ResidueTypes in the set have a "name3" tat matches the input name3
	bool has_name3( std::string const & name3 ) const;

	/// @brief Query a variant ResidueType by its base ResidueType and VariantType
	ResidueType const & get_residue_type_with_variant_added(
		ResidueType const & init_rsd,
		VariantType const new_type ) const;

	/// @brief return the residuetype we get from variant rsd type after removing the desired variant type
	ResidueType const & get_residue_type_with_variant_removed(
		ResidueType const & init_rsd,
		VariantType const old_type ) const;

	/// @brief query ResidueTypes by their AA enum type
	///
	/// @details similar to name3_map, return all matched residue types
	/// or an empty list.
	ResidueTypeCOPs const &
	aa_map_DO_NOT_USE( AA const & aa ) const;

	/// @brief select a set of ResidueTypes give certain criteria
	void
	select_residues_DO_NOT_USE(
		ResidueTypeSelector const & selector,
		ResidueTypeCOPs & matches
	) const;

	/// @brief beginning of aas_defined_ list
	std::list< AA >::const_iterator
	aas_defined_begin() const;

	/// @brief end of aas_defined_ list
	std::list< AA >::const_iterator
	aas_defined_end() const;

	const_residue_iterator
	all_residues_begin_DO_NOT_USE() const
	{
		return name_map_.begin();
	}

	const_residue_iterator
	all_residues_end_DO_NOT_USE() const
	{
		return name_map_.end();
	}

	/// alternate access to all residuetypes as vector
	ResidueTypeCOPs const &
	residue_types_DO_NOT_USE() const
	{
		return residue_types_;
	}

	/// @brief accessor for database_directory
	std::string const&
	database_directory() const
	{
		return database_directory_;
	}

	/// @brief the residues with no patches applied
	ResidueTypeCOPs base_residue_types() const { return base_residue_types_; }

	/// @brief the patches
	utility::vector1< PatchCOP > const & patches() const { return patches_; }

	/// @brief the patches, index by name.
	std::map< std::string, utility::vector1< PatchCOP > > const & patch_map() const { return patch_map_; }

	/// @brief apply patch -- look in ResidueTypeSet if it exists already.
	ResidueTypeCOP
	get_residue_type_with_patch( PatchCOP patch, ResidueTypeCOP rsd_type ) const;

	//////////////////
	// private methods
private:

	/// @brief clear residue maps
	void
	clear_residue_maps();

	/// @brief update residue maps
	void
	update_residue_maps();

	void
	add_residue_type_to_maps( ResidueTypeCOP rsd );

	void
	remove_residue_type_from_maps(ResidueTypeCOP rsd);

	bool
	make_sure_instantiated( ResidueTypeCOP const & rsd_type ) const;

	bool
	make_sure_instantiated( ResidueTypeCOPs const & rsd_types ) const;

	void
	figure_out_last_patch_from_name( std::string const & rsd_name,
		std::string & rsd_name_base,
		std::string & patch_name ) const;

	/// @brief helper function used during replacing residue types after, e.g., orbitals.
	void
	update_base_residue_types_if_replaced( ResidueTypeCOP rsd_type, ResidueTypeCOP rsd_type_new );

	/// @brief replace residue type -- pointer will stay the same, but the residue type information will change; useful for instantiating placeholders.
	void
	replace_residue_type_in_set_defying_constness( ResidueTypeCOP rsd_type, ResidueType const & rsd_new ) const;

	//////////////////
	// data
private:

	/// What does the ChemicalManager call this ResidueTypeSet?
	std::string name_;

	// The default subsidiary typesets, typically specified in the database summary file.
	// You can add a residue type with a different subsidiary typeset, but you'll have to
	// construct it yourself.
	AtomTypeSetCOP atom_types_;
	ElementSetCOP elements_;
	MMAtomTypeSetCOP mm_atom_types_;
	orbitals::OrbitalTypeSetCOP orbital_types_;

	/// @brief the residues
	ResidueTypeCOPs residue_types_;

	/// @brief null list of residues when query fails
	//should make this static or something
	ResidueTypeCOPs empty_residue_list_;

	/// @brief map to ResidueType pointers by AA enum
	std::map< AA, ResidueTypeCOPs > aa_map_;

	/// @brief map between ResidueType's interchangeability group string and a vector of ResidueTypeCOPs
	std::map< std::string, ResidueTypeCOPs > interchangeability_group_map_;

	/// @brief map to ResidueType pointers by 3-letter string name
	std::map< std::string, ResidueTypeCOPs > name3_map_;

	/// @brief map to ResidueType pointers by unique residue id
	std::map< std::string, ResidueTypeCOP > name_map_;

	/// @brief caching queries based on aa & variants to avoid recomputation with ResidueTypeFinder
	mutable std::map< std::pair< AA, utility::vector1< std::string > >, ResidueTypeCOPs > cached_aa_variants_map_;

	/// @brief map to ResidueType pointers that were replaced during application
	//  Do not delete yet -- may revive later to 'properly' handle replace_residue_types -- rhiju.
	// std::map< std::string, ResidueTypeCOP > replaced_name_map_;

	/// @brief list of AA types defined
	std::list< AA > aas_defined_;

	/// @brief the database directory of the generating files ---> allows to use cached dunbrack libs
	const std::string database_directory_;

	/// @brief the residues with no patches applied
	ResidueTypeCOPs base_residue_types_;
	std::map< std::string, std::set< std::string > > name3_generated_by_base_residue_name_;

	/// @brief the patches
	utility::vector1< PatchCOP > patches_;
	std::map< std::string, utility::vector1< PatchCOP > > patch_map_;

	bool on_the_fly_;

private:
	// uncopyable
	ResidueTypeSet( ResidueTypeSet const & );
	ResidueTypeSet const & operator = ( ResidueTypeSet const & );
};

} // chemical
} // core


#endif
