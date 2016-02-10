// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/Record.hh
/// @brief  Typedef for a .pdb file record data structure.
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// @author Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_io_pdb_Record_HH
#define INCLUDED_core_io_pdb_Record_HH

// Unit header
#include <core/io/pdb/Field.fwd.hh>

// C++ headers
#include <map>
#include <string>
#include <iostream>


namespace core {
namespace io {
namespace pdb {

/// @brief  A type for storing a .pdb file record as a map of field names to Fields.
typedef std::map< std::string, Field > Record;

// Insertion operator (overloaded so that Record can be "printed").
std::ostream & operator<<( std::ostream & os, Record const & record );

}  // namespace pdb
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_pdb_Record_HH
