// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/GlobalResidueTypeSet.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_GlobalResidueTypeSet_hh
#define INCLUDED_core_chemical_GlobalResidueTypeSet_hh


// Unit headers
#include <core/chemical/GlobalResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueTypeSet.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSetCache.fwd.hh>

#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Metapatch.fwd.hh>
#include <core/chemical/Patch.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>
#include <set>

namespace core {
namespace chemical {

/// @brief A collection of ResidueType defined

class GlobalResidueTypeSet : public ResidueTypeSet
{

public:

	/// @brief default c-tor
	GlobalResidueTypeSet() = delete; // Need to have directory for a Global RTS.

	/// @brief constructor from directory
	GlobalResidueTypeSet(
		std::string const & name,
		std::string const & directory
	);

	virtual ~GlobalResidueTypeSet();

	////////////////////////////
	// Private Initializer Functions
private:

	/// @brief Read in params files from the given Rosetta database directory.
	void init_restypes_from_database();

	/// @brief Read in the extra resdiue types specified with command line options.
	void init_restypes_from_commandline();

	/// @brief Provide a list of params files based on commandline options, appropriate for the current RTS mode
	utility::vector1< std::string > params_files_from_commandline() const;

	/// @brief Load various non-param files from commandline.
	/// These will be loaded as full atom types - it's the responsibility of the caller to convert.
	utility::vector1< ResidueTypeOP > extra_nonparam_restypes_from_commandline() const;

	/// @brief Load residue types from the command line-specified SQL database
	void load_residue_types_from_sql_database();

	/// @brief Read in patch and metapatch files from the given Rosetta database directory.
	void init_patches_from_database();

	/// @brief Read in the extra patches specified with command line options.
	void init_patches_from_commandline();

	/// @brief Deal with some special cases in patching
	/// @details Replacing residue types and D_AAs
	void deal_with_patch_special_cases();

	/// @brief apply patches to base ResidueType to generate variant ResidueTyes
	void
	place_adducts();

	//////////////////
	// public methods
public:

	/// @brief name of the residue type set (may be empty)
	/// @details The difference between a ResidueTypeSet *name* and a ResidueTypeSet *category* is that a
	/// a ResidueTypeSet *name* should uniquely identify a ResidueTypeSet (at lease those within the ChemicalManger)
	/// but more than one ResidueTypeSet may have the same *category*.
	/// The type specifies what compatibility class (full atom, centroid) the ResidueTypeSet has.
	/// Generally speaking, the *name* should only be used when interacting with the user.
	std::string const &
	name() const {
		return name_;
	}

	/// @brief query if a ResidueType of the unique residue id (name) is present.
	bool has_name( std::string const & name ) const override;

	/// @brief Gets all types with the given aa type and variants
	/// @details The number of variants must match exactly. Variants can be custom variants.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	ResidueTypeCOPs
	get_all_types_with_variants_aa( AA aa, utility::vector1< std::string > const & variants ) const override;

	/// @brief Gets all types with the given aa type and variants, making exceptions for some variants.
	/// @details The number of variants must match exactly. Variants can be custom variants, but exceptions must
	///           be standard types, listed in VariantType.hh.
	/// (It's assumed that the passed VariantTypeList contains no duplicates.)
	ResidueTypeCOPs
	get_all_types_with_variants_aa( AA aa,
		utility::vector1< std::string > const & variants,
		utility::vector1< VariantType > const & exceptions ) const override;

	/// @brief accessor for database_directory
	std::string const&
	database_directory() const
	{
		return database_directory_;
	}

	//////////////////
	// protected methods, implement from base class
protected:

	bool
	has_name_write_locked( std::string const & name ) const override;

	/// @brief Attempt to lazily load the given residue type from data.
	bool
	lazy_load_base_type( std::string const & rsd_base_name ) const override;

	/// @brief delete an unpatchable residue type from the set.
	/// This should never be called in a GlobalResidueTypeSet - it's here just to match the base class interface.
	virtual
	void
	remove_base_residue_type( std::string const & name ) override;

	/// @brief delete an unpatchable residue type from the set.
	/// This should never be called in a GlobalResidueTypeSet - it's here just to match the base class interface.
	virtual
	void
	remove_unpatchable_residue_type( std::string const & name ) override;

	//////////////////
	// private methods
private:

	// Set the name of the ResidueTypeSet
	void
	name( std::string const & setting ) { name_ = setting; }

	void
	pdb_components_filename( std::string const & setting ) { pdb_components_filename_ = setting; }

	std::string const &
	pdb_components_filename() const { return pdb_components_filename_; }

	void
	generate_all_residue_types();

	/// @brief From a file, read which IDs shouldn't be loaded from the components.
	void
	load_shadowed_ids( std::string const & directory, std::string const & file = "shadow_list.txt" );

	/// @brief Load a residue type from the components dictionary.
	ResidueTypeOP
	load_pdb_component( std::string const & pdb_id ) const;

	//////////////////
	// data
private:

	/// What does the ChemicalManager call this GlobalResidueTypeSet?
	std::string name_;

	/// @brief the database directory of the generating files ---> allows to use cached dunbrack libs
	const std::string database_directory_;

	/// @brief Which components shouldn't be loaded from the components file.
	std::set< std::string > shadowed_ids_;

	/// @brief data for lazy loading of PDB components
	std::string pdb_components_filename_;

private:
	// uncopyable
	GlobalResidueTypeSet( GlobalResidueTypeSet const & ) = delete;
	GlobalResidueTypeSet const & operator = ( GlobalResidueTypeSet const & ) = delete;

	// Regarding SERIALIZATION, the deleted default constructor should keep us from ever
	// accidentally serializing a GlobalResidueTypeSet.
};

} // chemical
} // core


#endif
