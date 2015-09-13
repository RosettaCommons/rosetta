// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/carbohydrates/CarbohydrateInfoManager.hh
/// @brief   Method declarations for CarbohydrateInfoManager.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_carbohydrates_CarbohydrateInfoManager_HH
#define INCLUDED_core_chemical_carbohydrates_CarbohydrateInfoManager_HH

// Unit header
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/SingletonBase.hh>

// C++ header
#include <string>
#include <map>


namespace core {
namespace chemical {
namespace carbohydrates {

/// @details  This class is a singleton and manages CarbohydratesInfo data that should only be read from the database
/// one time and shared among all instances of CarbohydrateInfo.
class CarbohydrateInfoManager : public utility::SingletonBase< CarbohydrateInfoManager > {
public:  // Declare friends ///////////////////////////////////////////////////
	friend class utility::SingletonBase< CarbohydrateInfoManager >;


public:  // Static constant data access ///////////////////////////////////////
	/// @brief Is the given 3-letter code a valid Rosetta/IUPAC code for carbohydrates?
	static bool is_valid_sugar_code( std::string const & code );

	/// @brief Get the monosaccharide root name from the given Rosetta/IUPAC 3-letter code.
	static std::string const & root_from_code( std::string const & code );

	/// @brief Get the 1-letter affix for designating a carbohydrate ring of this size.
	static char ring_affix_from_ring_size( core::Size ring_size );

	/// @brief Get the morpheme for designating a carbohydrate ring of this size.
	static std::string const & morpheme_from_ring_size( core::Size ring_size );


private:  // Private methods //////////////////////////////////////////////////
	// Empty constructor
	CarbohydrateInfoManager();

	// Singleton-creation function for use with utility::thread::threadsafe_singleton
	static CarbohydrateInfoManager * create_singleton_instance();

	// Get the map of Rosetta PDB 3-letter codes for saccharide residues mapped to the corresponding root requested,
	// creating it if necessary.
	// Called by the public static method root_from_code().
	std::map< std::string, std::string > const & code_to_root_map();

	// Get the map of carbohydrate ring sizes and their 1-letter affixes and morphemes requested, creating it if
	// necessary.
	// Called by the public static methods ring_affix_from_ring_size() and morpheme_from_ring_size().
	std::map< core::Size, std::pair< char, std::string > > const & ring_size_to_morphemes_map();


private:  // Private data /////////////////////////////////////////////////////
	std::map< std::string, std::string > code_to_root_map_;
	std::map< core::Size, std::pair< char, std::string > > ring_size_to_morphemes_map_;
};

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_carbohydrates_CarbohydrateInfoManager_HH
