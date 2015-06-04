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

#ifndef INCLUDED_core_io_pdb_Field_hh
#define INCLUDED_core_io_pdb_Field_hh


// Unit headers
#include <core/io/pdb/Field.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/Show.hh>

// c++ headers
#include <iostream>
#include <string>

namespace core {
namespace io {
namespace pdb {


/// @brief Data type Class to represent one field in PDB file.
class Field : public utility::Show {
public:

	Field();
	Field(Size s, Size e);
	Field(std::string type_, Size s, Size e);

	/// @brief read field value from given string.
	void getValueFrom(std::string source);

/// This class is intended to be just 'data' type class
/// no need to make it private.
public:

	/// @brief string value of field, type of the field.
	std::string type, value;

	/// @brief beginning position in line, ending position in line
	Size start, end;

	/// @brief Debug output.
	friend
	std::ostream&
	operator <<(std::ostream &os, Field const & F);

	/// @brief collection builder
	static RecordRef & getRecordCollection();

private:

};


/// @brief Debug printing, serialazing to Tracer like object.
std::ostream&
operator <<(std::ostream &os, Record const & R);

} // pdb
} // io
} // core

#endif // INCLUDED_core_io_pdb_Field_hh
