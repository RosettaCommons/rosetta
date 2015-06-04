// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/RingConformerManager.hh
/// @brief   Declarations and simple accessor/mutator definitions for NomenclatureManager.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_io_pdb_NomenclatureManager_HH
#define INCLUDED_core_io_pdb_NomenclatureManager_HH

// Unit header
#include <core/io/pdb/NomenclatureManager.fwd.hh>

// Package header
#include <core/io/pdb/alt_codes_io.hh>

// Utility header
#include <utility/SingletonBase.hh>

// C++ header
#include <string>


namespace core {
namespace io {
namespace pdb {

/// @details  This class is a singleton and manages AltCodeMap data that should only be read from the database one time
/// and shared among all processes constructing Poses.
class NomenclatureManager : public utility::SingletonBase< NomenclatureManager > {
public:  // Declare friends ///////////////////////////////////////////////////
	friend class utility::SingletonBase< NomenclatureManager >;


public:  // Static constant data access ///////////////////////////////////////
	/// @brief  Return a pair of Rosetta names (3-letter code and base ResidueType name, if available) from the
	/// given PDB 3-letter code.
	static std::pair< std::string, std::string > rosetta_names_from_pdb_code( std::string const & pdb_code );


private:  // Private methods //////////////////////////////////////////////////
	// Empty constructor
	NomenclatureManager();

	// Singleton-creation function for use with utility::thread::threadsafe_singleton
	static NomenclatureManager * create_singleton_instance();

	// Get the map requested, creating it if necessary.
	// Called by the public static method rosetta_names_from_pdb_code()
	AltCodeMap const & get_alternate_3_letter_code_map();

	// Try various combinations to locate the specific file being requested by the user.
	// (inspired by core::scoring::ScoreFunction::find_weights_file())
	std::string find_alternate_codes_file( std::string filename );

private:  // Private data /////////////////////////////////////////////////////
	AltCodeMap alt_codes_;
};

}  // namespace pdb
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_pdb_NomenclatureManager_HH

