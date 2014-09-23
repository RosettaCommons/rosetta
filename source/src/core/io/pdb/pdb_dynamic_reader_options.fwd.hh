// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/pdb_dynamic_reader_options.fwd.hh
///
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_core_io_pdb_pdb_dynamic_reader_options_FWD_HH
#define INCLUDED_core_io_pdb_pdb_dynamic_reader_options_FWD_HH


#include <utility/pointer/owning_ptr.hh>

// C++ headers

namespace core {
namespace io {
namespace pdb {

class PDB_DReaderOptions;

typedef utility::pointer::shared_ptr< PDB_DReaderOptions > PDB_DReaderOptionsOP;
typedef utility::pointer::shared_ptr< PDB_DReaderOptions const > PDB_DReaderOptionsCOP;

} // namespace pdb
} // namespace io
} // namespace core


#endif // INCLUDED_core_io_pdb_pdb_dynamic_reader_options_FWD_HH
