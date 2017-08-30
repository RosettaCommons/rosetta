// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/carbohydrates/CarbohydrateInfoManager.hh
/// @brief   Method declarations for CarbohydrateInfoManager.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @author  Vikram K. Mulligan (vmullig@uw.edu) -- Made the CarbohydrateInfoManager threadsafe.


#ifndef INCLUDED_core_chemical_carbohydrates_CarbohydrateInfoManager_HH
#define INCLUDED_core_chemical_carbohydrates_CarbohydrateInfoManager_HH

// Unit header
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.fwd.hh>

// Package header
#include <core/chemical/carbohydrates/SugarModificationsNomenclatureTable.hh>
#include <core/chemical/carbohydrates/LinkageConformers.hh>

// Project headers
#include <core/types.hh>

#include <core/chemical/VariantType.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>
#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif

// C++ header
#include <string>
#include <map>


namespace core {
namespace chemical {
namespace carbohydrates {

struct RootData {
	std::string root; // The root of the name
	char DL; // The D/L specification
	char anomeric; // The location of the anomeric carbon (aldose/2-ketose/3-ketose, etc.)
};

/// @details  This class is a singleton and manages CarbohydrateInfo data that should only be read from the database
/// one time and shared among all instances of CarbohydrateInfo.
class CarbohydrateInfoManager : public utility::SingletonBase< CarbohydrateInfoManager > {
public:  // Declare friends ///////////////////////////////////////////////////
	friend class utility::SingletonBase< CarbohydrateInfoManager >;


public:  // Static constant data access ///////////////////////////////////////
	/// @brief  Is the given 3-letter code a valid Rosetta/IUPAC code for carbohydrates?
	static bool is_valid_sugar_code( std::string const & code );

	/// @brief  Get the monosaccharide root name from the given Rosetta/IUPAC 3-letter code.
	static std::string const & root_from_code( std::string const & code );

	/// @brief  Get the default stereochemistry from the given Rosetta/IUPAC 3-letter code.
	static char default_stereochem_from_code( std::string const & code );

	/// @brief  Get the position of the anomeric carbon from the given Rosetta/IUPAC 3-letter code.
	static char anomeric_position_from_code( std::string const & code );

	/// @brief  Is the given 1-letter code valid for designating carbohydrate ring size?
	static bool is_valid_ring_affix( char affix );

	/// @brief  Get the 1-letter affix for designating a carbohydrate ring of this size.
	static char ring_affix_from_ring_size( core::Size ring_size );

	/// @brief  Get the morpheme for designating a carbohydrate ring of this size.
	static std::string const & morpheme_from_ring_size( core::Size ring_size );


	/// @brief  Is the given short affix valid for designating a sugar modification?
	static bool is_valid_modification_affix( std::string const & affix );

	/// @brief  Get the Rosetta patch name for this sugar modification affix.
	static std::string const & patch_name_from_affix( std::string const & affix );

	/// @brief  Get the default position for this sugar modification affix.
	static core::uint default_position_from_affix( std::string const & affix );

	/// @brief  Does this sugar modification affix have an inherent/defined position?
	static bool affix_has_inherent_position( std::string const & affix );


	/// @brief  Get the Cn_BRANCH_POINT VariantType for this atom name, e.g., On.
	static VariantType branch_variant_type_from_atom_name( std::string const & atom_name );

	/// @brief  Get the Cn_BRANCH_POINT VariantType for this position (n).
	static VariantType branch_variant_type_from_position( core::uint const position );


	/// @brief  Does the linkage between the given pair of monosaccharide residues have statistics in the database?
	static bool pair_has_linkage_statistics( std::string const & res1, std::string const & res2 );

	/// @brief  Get the "linkage conformer" statistical data from a given pair of monosaccharide residues.
	static utility::vector1< LinkageConformerData > linkages_from_pair(
		std::string const & res1, std::string const & res2 );

	/// @brief Get a map of short names to the full iupac glycan sequence for common glycosylations.
	static std::map< std::string, std::string> const & get_short_name_to_iupac_strings_map();

private:  // Private methods //////////////////////////////////////////////////
	// Empty constructor
	CarbohydrateInfoManager();

	// Get the map of Rosetta PDB 3-letter codes for saccharide residues mapped to the corresponding root and default
	// stereochemistry, creating it if necessary.
	// Called by the public static method root_from_code() and default_stereochem_from_code().
	std::map< std::string, RootData > const & code_to_root_map();


	// Get the map of carbohydrate ring sizes and their 1-letter affixes and morphemes, creating it if
	// necessary.
	// Called by the public static methods ring_affix_from_ring_size() and morpheme_from_ring_size().
	std::map< core::Size, std::pair< char, std::string > > const & ring_size_to_morphemes_map();

	// Get a list of valid 1-letter affixes for ring size, creating it if necessary.
	utility::vector1< char > const & ring_affixes();


	// Get the table of nomenclature data for sugar modifications, creating it if necessary.
	SugarModificationsNomenclatureTable const & nomenclature_table();

	// Get a map of sugar modification affixes to Rosetta patch names, creating it if necessary.
	std::map< std::string, std::string > const & affix_to_patch_map();

	// Get a map of sugar modification affixes to default positions, creating it if necessary.
	std::map< std::string, core::uint > const & affix_to_position_map();

	// Get a map of sugar modification affixes to a boolean indication of if the position is inherent,
	// creating it if necessary.
	std::map< std::string, bool > const & affix_to_position_inherency_map();


	// Get a map of linkage conformer statistical data, creating it if necessary.
	LinkageConformers const & linkage_conformers_map();


	// Try various combinations to locate the specific file being requested by the user.
	// (inspired by core::scoring::ScoreFunction::find_weights_file())
	std::string find_linkage_conformer_data_file( std::string filename );

	/// @brief Get a map of short names to the full iupac glycan sequence for common glycosylations.
	std::map< std::string, std::string> const & short_name_to_iupac_strings_map();

private:  // Private data /////////////////////////////////////////////////////

#ifdef MULTI_THREADED
	utility::thread::ReadWriteMutex code_to_root_map_mutex_;
	utility::thread::ReadWriteMutex ring_size_to_morphemes_mutex_;
	utility::thread::ReadWriteMutex ring_affixes_mutex_;
	utility::thread::ReadWriteMutex nomenclature_table_mutex_;
	utility::thread::ReadWriteMutex affix_to_patch_mutex_;
	utility::thread::ReadWriteMutex affix_to_position_mutex_;
	utility::thread::ReadWriteMutex affix_to_position_inherency_mutex_;
	utility::thread::ReadWriteMutex linkage_conformers_mutex_;
	utility::thread::ReadWriteMutex short_name_to_iupac_strings_mutex_;
#endif

	std::map< std::string, RootData > code_to_root_map_;  // also stores default stereochemistry

	std::map< core::Size, std::pair< char, std::string > > ring_size_to_morphemes_map_;
	utility::vector1< char > ring_affixes_;

	// Nomenclature
	SugarModificationsNomenclatureTable nomenclature_table_;
	std::map< std::string, std::string > affix_to_patch_map_;
	std::map< std::string, core::uint > affix_to_position_map_;
	std::map< std::string, std::string > short_name_to_iupac_strings_map_;
	std::map< std::string, bool > affix_to_position_inherency_map_;

	// Glycan Relax
	std::map< std::pair< std::string, std::string >, utility::vector1< LinkageConformerData > >
		linkage_conformers_map_;
};  // class CarbohydrateInfoManager


// Helper function ////////////////////////////////////////////////////////////
std::pair< std::string, std::string > convert_residue_names_into_linkage_map_key(
	std::string const & name1, std::string const & name2 );

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_carbohydrates_CarbohydrateInfoManager_HH
