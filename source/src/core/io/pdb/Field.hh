// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/Field.hh
/// @brief  Definitions for the Field data structures and related helper function declarations.
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// @author Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_io_pdb_Field_HH
#define INCLUDED_core_io_pdb_Field_HH

// Unit header
#include <core/io/pdb/Field.fwd.hh>

// Project header
#include <core/types.hh>

// C++ headers
#include <iostream>
#include <string>


namespace core {
namespace io {
namespace pdb {

struct Field {
	// Methods ////////////////////////////////////////////////////////////////
	// Empty constructor
	Field();

	/// @brief  Convenience constructor.
	Field( core::uint start_in, core::uint end_in );


	/// @brief  Read field value from given .pdb line and set.
	void set_value_from_pdb_line( std::string source );


	// Member data ////////////////////////////////////////////////////////////
	/// @brief  String value of field.
	std::string value;

	/// @brief  Beginning position in line, ending position in line.
	core::uint start, end;
};  // struct Field


// Helper function ////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that Field can be "printed").
std::ostream & operator<<( std::ostream & os, Field const & field );

}  // pdb
}  // io
}  // core

#endif  // INCLUDED_core_io_pdb_Field_HH
