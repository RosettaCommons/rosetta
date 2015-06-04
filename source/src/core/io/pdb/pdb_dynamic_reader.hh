// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/pdb_dynamic_reader.hh
/// @brief  Declarations for PDB Dynamic reader.
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)

// Note: DO NOT ACCESS THE OPTIONS SYSTEM DIRECTLY IN THIS FILE!
// Doing so will mean the Resource Manager will not work properly.
// Instead, modify PDB_DReaderOptions to include the option


#ifndef INCLUDED_core_io_pdb_pdb_dynamic_reader_hh
#define INCLUDED_core_io_pdb_pdb_dynamic_reader_hh


// Unit headers
#include <core/io/pdb/file_data.fwd.hh>
#include <core/io/pdb/pdb_dynamic_reader_options.fwd.hh>

// Project headers
#include <core/io/pdb/Field.hh>

// Utility headers
#include <utility/Show.hh>
#include <utility/vector1.fwd.hh>

// C++ headers
#ifdef WIN32
#include <string>
#endif

#include <map>
#include <vector>
#include <iostream>


namespace core {
namespace io {
namespace pdb {


typedef std::string String;

/// @brief PDB Reader itself, D - for dynamic approach of type handling
class PDB_DReader
{
public:
	/// @brief creating record from given string. Also, read Field values from string.
	static Record mapStringToRecord(const String & s);


	/// @brief Reverse operation - create PDB string from given Record
	static String createPDBString(const Record & R);


	/// @brief Parse whole PDB string and return vector of records in order they was in PDB.
	static std::vector<Record> parse(const String &);


	/// @brief create File data structure from array of Records.
	static FileData createFileData(std::vector<Record> &);

	/// @brief create File data structure from array of Records and a set of options.
	static FileData createFileData(std::vector<Record> &, PDB_DReaderOptions const & options);

	/// @brief create File data structure from string containing PDB information.
	static FileData createFileData(const String & data);

	/// @brief create File data structure from string containing PDB information and a set of options.
	static FileData createFileData(const String & data, PDB_DReaderOptions const & options);
	
	/// @brief create PDB-like string to represent given FileData object
	static String createPDBData(FileData const &fd);

	/// @brief create PDB-like vector of strings to represent given FileData object.
	static utility::vector1< std::string > createPDBData_vector(FileData const & fd );


	/// @brief create vector of records for given FileData object.
	static std::vector<Record> createRecords(FileData const & fd);

};

/// @brief print int with format to string
std::string print_i(const char *format, int I);

/// @brief print double with format to string
std::string print_d(const char *format, double d);

} // namespace pdb
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_pdb_pdb_dynamic_reader_HH
