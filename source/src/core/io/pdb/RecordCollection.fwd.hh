// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/pdb/RecordCollection.fwd.hh
/// @brief   Forward declaration for RecordCollection.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_io_pdb_RecordCollection_FWD_HH
#define INCLUDED_core_io_pdb_RecordCollection_FWD_HH

// Unit headers
#include <core/io/pdb/Record.hh>
#include <core/io/pdb/RecordType.hh>

// C++ header
#include <map>


namespace core {
namespace io {
namespace pdb {

/// @brief  A singleton class for handling static const PDB Record definition data.
class RecordCollection;

/// @brief  Definitions of all possible records (line types), that can exist in a PDB file.
typedef std::map< RecordType, Record > RecordDef;

}  // namespace pdb
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_pdb_RecordCollection_FWD_HH
