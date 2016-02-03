// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/Field.hh
/// @brief Each line of a PDB file is a Record which is divided into Fields
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// @author Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_io_pdb_Field_HH
#define INCLUDED_core_io_pdb_Field_HH

// Unit headers
#include <core/io/pdb/Field.fwd.hh>
#include <core/io/pdb/PDBDataType.hh>

// Project header
#include <core/types.hh>

// C++ headers
#include <iostream>
#include <string>


namespace core {
namespace io {
namespace pdb {

/// @brief Data type Class to represent one field in PDB file.
class Field {
public:

	Field();
	Field( core::uint start_in, core::uint end_in, PDBDataType data_type_in );

	/// @brief Read field value from given string and set.
	void set_value_from_string( std::string source );

	/// This class is intended to be just 'data' type class
	/// no need to make it private.
public:
	/// @brief String value of field
	std::string value;

	/// @brief PDB-defined data type for this field
	PDBDataType data_type;

	/// @brief beginning position in line, ending position in line
	core::uint start, end;

	/// @brief collection builder
	//static RecordRef & getRecordCollection();
};

// Helper functions
/// @brief Get the PDBDataType value from the corresponding string.
PDBDataType get_pdb_data_type_from_string( std::string const & type );

/// @brief Get the string from the corresponding PDBDataType value.
std::string get_pdb_data_type_from_string( PDBDataType type );

/// @brief Debug output.
std::ostream & operator<<(std::ostream & os, Field const & F );

/// @brief Debug printing, serializing to Tracer like object.
std::ostream & operator<<( std::ostream & os, Record const & R );

}  // pdb
}  // io
}  // core

#endif // INCLUDED_core_io_pdb_Field_HH
