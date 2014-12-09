// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/Field.fwd.hh
/// @brief Each line of a PDB file is a Record which is divided into Fields
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_core_io_pdb_Field_fwd_hh
#define INCLUDED_core_io_pdb_Field_fwd_hh

#include <map>
#include <string>

namespace core {
namespace io {
namespace pdb {


// forward declaration
class Field;

typedef std::map<std::string, Field> Record;

/// @brief collection of all possible records (line types), that can exist in a PDB file.
typedef std::map<std::string, Record> RecordRef;


} // pdb
} // io
} // core

#endif // INCLUDED_core_io_pdb_Field_fwd_hh
