// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/carbohydrates/database_io.hh
/// @brief   Database input/output function declarations for carbohydrate-specific data.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_carbohydrates_database_io_HH
#define INCLUDED_core_chemical_carbohydrates_database_io_HH

// Package header
#include <core/chemical/carbohydrates/SugarModificationsNomenclatureTable.hh>
#include <core/chemical/carbohydrates/carbohydrate_data_structures.hh>

// Project header
#include <core/types.hh>

// Utility header
#include <utility/vector1.hh>

// C++ headers
#include <map>
#include <iosfwd>


namespace core {
namespace chemical {
namespace carbohydrates {

/// @brief  Try various combinations to locate the specific glycan sequence file being requested by the user.
std::string find_glycan_sequence_file( std::string filename );

/// @brief  Read a single-line glycan sequence file.
std::string read_glycan_sequence_file( std::string filename );


/// @brief  Return a map of strings to strings, which are saccharide-specific 3-letter codes mapped to IUPAC roots, read
/// from a database file.
std::map< std::string, std::string >
read_codes_and_roots_from_database_file( std::string const & filename );

/// @brief  Return a map of Sizes to pairs of char and string, which are ring sizes mapped to 1-letter affixes and
/// morphemes, respectively, read from a database file.
std::map< core::Size, std::pair< char, std::string > >
read_ring_sizes_and_morphemes_from_database_file( std::string const & filename );

/// @brief  Return a table of nomenclature data for sugar modifications, read from a database file.
SugarModificationsNomenclatureTable
read_nomenclature_table_from_database_file( std::string const & filename );

/// @brief  Return a map of linkage conformer data, read from a database file.
LinkageConformers
read_linkage_conformers_from_database_file( std::string const & filename );

///@breif  Return a map of short names to IUPAC formatted strings.
///  Reads from db_dir/common_names.txt, and loads the IUPAC files as strings.
std::map< std::string, std::string >
read_short_names_to_iupac_format_string( std::string const & dir, std::string common_mapping_path);



}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_carbohydrates_database_io_HH
