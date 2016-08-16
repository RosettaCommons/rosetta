// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/io/pdb/pdb_reader.hh
/// @brief   Function declarations for reading of .pdb files.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @note    We may perhaps decide in the future to wrap this functionality in
/// a class.  ~Labonte


#ifndef INCLUDED_core_io_pdb_pdb_reader_HH
#define INCLUDED_core_io_pdb_pdb_reader_HH

// Unit headers
#include <core/io/pdb/Record.hh>

// Package headers
#include <core/io/StructFileRep.fwd.hh>
#include <core/io/StructFileReaderOptions.fwd.hh>  // TODO: Rename after refactoring is complete.

// Utility headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <string>


namespace core {
namespace io {
namespace pdb {

/// @brief  Convert a .pdb file line into a Record data structure.
Record create_record_from_pdb_line( std::string const & line );


/// @brief  Create a list of .pdb format records from the lines from a .pdb file.
utility::vector1< Record > create_records_from_pdb_lines( utility::vector1< std::string > const & lines );

/// @brief  Create a list of .pdb format records from the entire contents of a .pdb file.
utility::vector1< Record > create_records_from_pdb_file_contents( std::string const & pdb_contents );


/// @brief  Create a representation of structural file data from a list of .pdb format records with options.
StructFileRep create_sfr_from_pdb_records( utility::vector1< Record > & records, StructFileReaderOptions const & options );

/// @brief  Create a representation of structural file data from a list of .pdb format records.
StructFileRep create_sfr_from_pdb_records( utility::vector1< Record > & records );


/// @brief  Create a representation of structural file data from .pdb file contents with options.
StructFileRep create_sfr_from_pdb_file_contents( std::string const & pdb_contents, StructFileReaderOptions const & options );

/// @brief  Create a representation of structural file data from .pdb file contents.
StructFileRep create_sfr_from_pdb_file_contents( std::string const & pdb_contents );


/// @brief  Create a representation of structural file data from a .pdb file by file.
//StructFileRep create_sfr_from_pdb_file( utility::io::izstream const & file );

/// @brief  Create a representation of structural file data from a .pdb file by filename.
//StructFileRep create_sfr_from_pdb_file( std::string const & filename );


// .pdb Record Storage Functions //////////////////////////////////////////////
/// @brief  Convert .pdb SEQRES record into SFR data.
void store_chain_sequence_record_in_sfr( Record seqres_record, StructFileRep & sfr );

/// @brief  Convert .pdb MODRES record into SFR data.
void store_mod_res_record_in_sfr( Record modres_record, StructFileRep & sfr );


/// @brief  Parse .pdb HETNAM text field to extract full resID and convert into SFR data.
void store_base_residue_type_name_in_sfr(
	std::string const & hetID,
	std::string const & text_field,
	StructFileRep & sfr );

/// @brief  Convert .pdb HETNAM record into SFR data.
void store_heterogen_name_record_in_sfr( Record hetnam_record, StructFileRep & sfr );

/// @brief  Convert .pdb HETSYN record into SFR data.
void store_heterogen_synonym_record_in_sfr( Record hetsyn_record, StructFileRep & sfr );

/// @brief  Convert .pdb FORMUL record into SFR data.
void store_formula_record_in_sfr( Record formul_record, StructFileRep & sfr );


/// @brief  Convert .pdb SSBOND record into SFR data.
void store_ssbond_record_in_sfr( Record ssbond_record, StructFileRep & sfr );

/// @brief  Convert .pdb LINK record into SFR data.
void store_link_record_in_sfr( Record link_record, StructFileRep & sfr );

/// @brief  Convert .pdb CISPEP record into SFR data.
void store_cis_peptide_record_in_sfr( Record cispep_record, StructFileRep & sfr );


/// @brief  Convert .pdb CRYST1 record into SFR data.
void store_crystallographic_parameter_record_in_sfr( Record crystal_record, StructFileRep & sfr );


/// @brief  Parse and store unknown record types into SFR data.
void store_unknown_records_in_sfr( utility::vector1< Record > unknown_records, StructFileRep & sfr );

}  // namespace pdb
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_pdb_pdb_reader_HH
