// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/NomenclatureManager.hh
/// @brief   Declarations and simple accessor/mutator definitions for NomenclatureManager.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_io_NomenclatureManager_HH
#define INCLUDED_core_io_NomenclatureManager_HH

// Unit header
#include <core/io/NomenclatureManager.fwd.hh>

// Package header
#include <core/io/alt_codes_io.hh>

// Utility header
#include <utility/SingletonBase.hh>

// C++ header
#include <string>
#include <set>


namespace core {
namespace io {

/// @details  This class is a singleton and manages AltCodeMap data that should only be read from the database one time
/// and shared among all processes constructing Poses.
class NomenclatureManager : public utility::SingletonBase< NomenclatureManager > {
public:  // Declare friends ///////////////////////////////////////////////////
	friend class utility::SingletonBase< NomenclatureManager >;


public:  // Static constant data access ///////////////////////////////////////

	/// @brief  Return a pair of Rosetta names (3-letter code and base ResidueType name, if available) from the
	/// given PDB 3-letter code.
	static std::pair< std::string, std::string >
	rosetta_names_from_pdb_code( std::string const & pdb_code );

	static bool is_NA( std::string const & name3 );
	static bool is_old_RNA( std::string const & name3 );
	static bool is_old_DNA( std::string const & name3 );
	static bool decide_is_d_aa( std::string const & name3 );
	static bool decide_is_l_aa( std::string const & name3 );
	static bool decide_is_known_achiral( std::string const & name3 );
	static bool is_metal( std::string const & name3 );
	static bool is_sugar( std::string const & name3 );

	static std::string annotated_sequence_from_modomics_oneletter_sequence( std::string const & seq );
	static std::string annotated_sequence_from_IUPAC_sequence( std::string const & seq );

private:  // Private methods //////////////////////////////////////////////////
	// Empty constructor
	NomenclatureManager();

	// Get the map requested, creating it if necessary.
	// Called by the public static method rosetta_names_from_pdb_code()
	AltCodeMap const & get_alternate_3_letter_code_map() const;

	// Get the set requested, creating it if necessary.
	// Called by the public static method rosetta_names_from_pdb_code()
	std::set< std::string > const & na_set() const { return is_NA_; };
	std::set< std::string > const & old_rna_set() const { return is_old_RNA_; };
	std::set< std::string > const & old_dna_set() const { return is_old_DNA_; };
	std::set< std::string > const & d_aa_set() const { return d_aa_set_; };
	std::set< std::string > const & l_aa_set() const { return l_aa_set_; };
	std::set< std::string > const & achiral_set() const { return achiral_set_; };
	std::set< std::string > const & metal_set() const { return metal_set_; };
	std::set< std::string > const & sugar_set() const { return sugar_set_; };
	std::map< std::string, std::string > const & iupac_map() const { return annotated_seq_from_IUPAC_map_; };
	std::map< char, std::string > const & modomics_map() const { return annotated_seq_from_modomics_map_; };

	// Try various combinations to locate the specific file being requested by the user.
	// (inspired by core::scoring::ScoreFunction::find_weights_file())
	std::string find_alternate_codes_file( std::string const & filename );

private:  // Private data /////////////////////////////////////////////////////

	AltCodeMap alt_codes_;
	std::set< std::string > is_NA_;
	std::set< std::string > is_old_RNA_;
	std::set< std::string > is_old_DNA_;
	std::set< std::string > d_aa_set_;
	std::set< std::string > l_aa_set_;
	std::set< std::string > achiral_set_;
	std::set< std::string > metal_set_;
	std::set< std::string > sugar_set_;

	std::map< char, std::string > annotated_seq_from_modomics_map_;
	std::map< std::string, std::string > annotated_seq_from_IUPAC_map_;
};

}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_pdb_NomenclatureManager_HH
