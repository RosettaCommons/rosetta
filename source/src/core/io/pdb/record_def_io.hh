// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/pdb/record_def_io.hh
/// @brief   Database input/output function declarations for PDB record definition data.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_io_pdb_record_def_io_HH
#define INCLUDED_core_io_pdb_record_def_io_HH

// Unit headers
#include <core/io/pdb/RecordType.hh>
#include <core/io/pdb/RecordCollection.fwd.hh>

// C++ headers
#include <string>


namespace core {
namespace io {
namespace pdb {

/// @brief  Return a mapping of PDB record types to record definitions.
RecordDef read_record_definitions_from_file(
	std::string const & filename,
	std::map< std::string, RecordType > const & record_type_map );

}  // namespace pdb
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_pdb_record_def_io_HH
